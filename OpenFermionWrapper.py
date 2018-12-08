from openfermion import *
from openfermionpsi4 import run_psi4
# import os
# import time
# import numpy
# import sys
# from scipy.linalg import eigh

class OpenFermionWrapper():
	''' 
	This class is designed to act as a wrapper for the OpenFermion package.

	This wrapper's focus is to create an interface between a client program
	and the OpenFermion syntax with an increased ease of use. 

	The intended use case for this wrapper is when a client program would
	like to define a molecule and create the hamiltonian, energy eigenstates, 
	or a quantum circuit written in QASM of the time evolution of the
	hamiltonian through the Trotter-Suzuki decomposition method. 

	Given the increased complexity around these use-cases, this wrapper is 
	created to keep in mind that client programs may want to restructure the
	created data. Hence, in contrast to normal programming design, there will
	be no private members of this class.

	Typical Progression:
		Set variables using setter functions
		load_molecule
		create_hamiltonians
		perform_transform
		generate_circuit

	Author: William A. Simon
	Date: 12/7/2018

	TODO:
		support BKSF, BKTree, Parity mappings
		get N lowest energy states

	 '''

    def __init__(self):
        ''' Initialize empty instances '''
        self.name = None
        self.geometery = None
        self.basis = None
        self.multiplicity = None
        self.charge = None
        self.description = ''
        self.filename = None

        self.plugin = "psi4" # other option is "pyscf"

        self.molecule = None

        self.active_space_start = None
        self.active_space_stop = None

        self.molecular_hamiltonian = None
        self.fermion_hamiltonian = None
        self.n_qubits = None

        self.ground_state_energy = None

        self.qubit_hamiltonian_jw = None
        self.qubit_hamiltonian_bk = None
        self.qubit_hamiltonian_bksf = None
        self.qubit_hamiltonian_bktree = None
        
        self.jw_circuit = None
        self.bk_circuit = None

    def set_name(self, name): self.name = name
    def set_geometry(self, geometery): self.geometery = geometery
    def set_basis(self, basis): self.basis = basis
    def set_multiplicity(self, multiplicity): self.multiplicity = multiplicity
    def set_charge(self, charge): self.charge = charge
    def set_description(self, description): self.description = description
    def set_active_space_start(self, orbital): self.active_space_start = orbital
    def set_active_space_stop(self, orbital): self.active_space_stop = orbital
   	def set_plugin(self, plugin):
   		if((plugin != 'pyscf') and (plugin != "pyscf")):
   			sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
   								Plugin not recognized\n''')
   		self.plugin = plugin

    def mapping_error(self, mapping):
        ''' A simple helper function to print and exit with an error related
                to an incorrect mapping request. '''
        sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
        				 Mapping: {} was not recognized'''.format(mapping))

	def set_parameters(self, geometry, basis, multiplicity, charge, 
   		description):
		# Use setter functions in case client has not set from outside the
		#	  wrapper
   		set_geometry(geometry)
   		set_basis(basis)
   		set_multiplicity(multiplicity)
   		set_charge(charge)
   		set_description(description)

   	def check_parameters(self):
   		if((self.geometry == None)     or (self.basis == None) or
   		   (self.multiplicity == None) or (self.charge == None)):
   			sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
   								Must define all molecular parameters\n''')

   	def load_molecule(self, geometry=None, basis=None, multiplicity=None, 
   		charge=None, description="", forceCalculation=False):

   		set_parameters(geometry, basis, multiplicity, charge, description)
   		check_parameters()

        # Generate and populate instance of MolecularData data structure
        if(self.description != ''):
            molecule = MolecularData(self.geometry, self.basis, 
                                     self.multiplicity, self.charge,
                                     self.description)
        else:
            molecule = MolecularData(self.geometry, self.basis, 
                                     self.multiplicity, self.charge)

        # Determine if integrals have been previously generated or not
        if(not(os.path.exists(molecule.filename + '.hdf5')) or 
              (forceCalculation)):

        	if(self.plugin == 'psi4'):
	        	# Run Psi4 calculation protocol
	            run_scf, run_mp2, run_cisd, run_ccsd, run_fci = 1, 1, 1, 1, 1

	            molecule = run_psi4(molecule, run_scf=run_scf, run_mp2=run_mp2,
	                                run_cisd=run_cisd, run_ccsd=run_ccsd,
	                                run_fci=run_fci)
	        else:
	        	# Run PySCF calculation protocol
            
            # Save molecule for future, so that regeneration is not required
            molecule.save()
            molecule.load()
        else:
            molecule.load()

        self.molecule = molecule

    def create_hamiltonians(self):
        ''' This function is designed to generate the fermion hamiltonian and
            qubit hamiltonians for both the transformations allowed in 
            OpenFermion. '''

	    # Try to load the molecule if the client has not yet loaded it 
	    if(not self.molecule):
	    	load_molecule()

	    if((self.active_space_start) and (self.active_space_stop)):
		    # Get the Hamiltonian in an active space.
		    self.molecular_hamiltonian = 
		    	self.molecule.get_molecular_hamiltonian(
			        occupied_indices=range(self.active_space_start),
			        active_indices=range(self.active_space_start, 
			                             self.active_space_stop))
		else:
			self.molecular_hamiltonian = self.molecule.get_molecular_hamiltonian()
	    
	    self.fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
	    self.n_qubits = count_qubits(self.fermion_hamiltonian)

	def set_ground_state_energy(self):
		if(not self.fermion_hamiltonian):
			create_hamiltonians()

		sparse_matrix = get_sparse_operator(self.fermion_hamiltonian)
		self.ground_state_energy = get_ground_state(sparse_matrix)[0]

	def perform_transform(self, mapping):
		''' 
		This function performs the approprate transformation to the 
		FermionOperator representing the current molecule. 

		The currently supported mappings are:
			Bravyi-Kitaev ("BK")
			Jordan-Wigner ("JW")
		'''

		# Attempt to catch user logic error
		if(not self.fermion_hamiltonian):
			create_hamiltonians()

		if(mapping == "BK"):
			self.qubit_hamiltonian_bk = bravyi_kitaev(self.fermion_hamiltonian)
		elif(mapping == "JW"):
			self.qubit_hamiltonian_jw = jordan_wigner(self.fermion_hamiltonian)
		else:
			mapping_error(mapping)

	def generate_circuit(self, mapping=None):
		''' 
		This function uses the built-in cability in OpenFermion to create
		a quantum circuit in QASM that represents the exponentiation of the 
		hamiltonian over time using the Trotter-Suzuki decomposition method.

		In keeping with the theme, checks are made to attempt to catch user
		logic error and still result in success. Note: If a mapping is
		not defined, recursive calls are made to use all implemented 
		transformations. 
		'''
		if(mapping == 'BK'):
			# Attempt to catch user logic error
			if(not self.qubit_hamiltonian_bk):
				perform_transform(mapping)
			self.bk_circuit = pauli_exp_to_qasm([self.qubit_hamiltonian_bk])
		elif(mapping == 'JW'):
			# Attempt to catch user logic error
			if(not self.qubit_hamiltonian_jw):
				perform_transform(mapping)
			self.jw_circuit = pauli_exp_to_qasm([self.qubit_hamiltonian_jw])
		else:
			# Attempt to catch user logic error
			generate_circuit('BK')
			generate_circuit('JW')



