from generateCircuit import GenerateCircuit
import random
import math

MIN_LENGTH = 0.1;
MAX_LENGTH = 1;

# minimum bond angle is pi-MIN_THETA
# maximum bond angle is always 180
MIN_THETA = -math.pi/3;
MAX_THETA = math.pi/3;

MIN_PHI = -math.pi;
MAX_PHI = math.pi;

def createGeometry3D(bondLength, atoms, uniform, linear):
	# This function is designed to created the intended geometry
	# for a chain of Hydrogen atoms. The geometry is a list of 
	# tuples wherein the first value is a string describing the
	# molecule (H for hydrogen) and the second value is the 
	# coordinates for the atom in 3 dimensions. 
	index = 1;
	geometry = list();	

	# Add first H atom at 0, 0, 0
	geometry.append(tuple(("H", (0,0,0))));

	#initialize the bondLength 
	if(bondLength == "random"):
		bondLength = random.uniform(MIN_LENGTH, MAX_LENGTH);

	positionX = 0;
	positionY = 0;
	positionZ = 0;
	prevTheta = 0;
	# continue to add atoms until desired amount
	while(index < atoms):
		# Set the position for the next atom 
		if(uniform):
			if(linear):
				# Linear and uniform distance, so X coord just gets 
				# 	the preset bond length
				positionX += bondLength;
			else:
				# Non linear and uniform distance, so bond length 
				# 	stays the same, but theta is generated randomly 
				#	between 0 and 2Pi. Add to previous X and Y coords
				theta = random.uniform(MIN_THETA, MAX_THETA);
				theta += prevTheta
				phi = random.uniform(MIN_PHI, MAX_PHI);
				positionX += bondLength*math.sin(theta)*math.cos(phi);
				positionY += bondLength*math.sin(theta)*math.sin(phi);
				positionZ += bondLength*math.cos(theta);
				prevTheta = theta;
		else:
			# Non uniform, so generate a new bondlength in the range
			bondLength = random.uniform(MIN_LENGTH, MAX_LENGTH);
			if(linear):
				# linear, so just add bond length to previous X coord
				positionX += bondLength;
			else:
				# non linear, so again, generate theta randomly and 
				# 	update both X and Y coords with new bond length
				theta = random.uniform(MIN_THETA, MAX_THETA);
				theta += prevTheta;
				phi = random.uniform(MIN_PHI, MAX_PHI);
				positionX += bondLength*math.sin(theta)*math.cos(phi);
				positionY += bondLength*math.sin(theta)*math.sin(phi);
				positionZ += bondLength*math.cos(theta);
				prevTheta = theta;

		# create atom tuple and push it to the list
		atom = tuple(("H", (positionX, positionY, positionZ)));
		geometry.append(atom);
		index += 1;


	return geometry;

def createGeometry2D(bondLength, atoms, uniform, linear):
	# This function is designed to created the intended geometry
	# for a chain of Hydrogen atoms. The geometry is a list of 
	# tuples wherein the first value is a string describing the
	# molecule (H for hydrogen) and the second value is the 
	# coordinates for the atom in 3 dimensions. 
	index = 1;
	geometry = list();	

	# Add first H atom at 0, 0, 0
	geometry.append(tuple(("H", (0,0,0))));

	#initialize the bondLength 
	if(bondLength == "random"):
		bondLength = random.uniform(MIN_LENGTH, MAX_LENGTH);

	positionX = 0;
	positionY = 0;
	prevTheta = 0;
	# continue to add atoms until desired amount
	while(index < atoms):
		# Set the position for the next atom 
		if(uniform):
			if(linear):
				# Linear and uniform distance, so X coord just gets 
				# 	the preset bond length
				positionX += bondLength;
			else:
				# Non linear and uniform distance, so bond length 
				# 	stays the same, but theta is generated randomly 
				#	between 0 and 2Pi. Add to previous X and Y coords
				theta = random.uniform(MIN_THETA, MAX_THETA);
				theta += prevTheta
				positionX += bondLength*math.cos(theta);
				positionY += bondLength*math.sin(theta);
				prevTheta = theta;
		else:
			# Non uniform, so generate a new bondlength in the range
			bondLength = random.uniform(MIN_LENGTH, MAX_LENGTH);
			if(linear):
				# linear, so just add bond length to previous X coord
				positionX += bondLength;
			else:
				# non linear, so again, generate theta randomly and 
				# 	update both X and Y coords with new bond length
				theta = random.uniform(MIN_THETA, MAX_THETA);
				theta += prevTheta;
				positionX += bondLength*math.cos(theta);
				positionY += bondLength*math.sin(theta);
				prevTheta = theta;

		# create atom tuple and push it to the list
		atom = tuple(("H", (positionX, positionY, 0)));
		geometry.append(atom);
		index += 1;


	return geometry;


def generateHChain(bondLength=.7414, mapping="BK", atoms=2, uniform=True, 
				   linear=True, threedimensions=True):

	H2 = GenerateCircuit();
	H2.set_name("HChain");
	description = str(bondLength) + str(atoms) + str(uniform) + str(linear);
	H2.set_description(description);

	if(not threedimensions):
		geometry = createGeometry2D(bondLength, atoms, uniform, linear);
	else:
		geometry = createGeometry3D(bondLength, atoms, uniform, linear);

	H2.set_geometry(geometry);
	H2.set_basis("sto-3g");
	H2.set_multiplicity(1);
	H2.set_charge(atoms-2);
	H2.load_molecule();
	H2.create_hamiltonians();  # need to create fermion hamiltonian
	H2.create_circuits(mapping);  # maps fermion to qubit hamil and circuit
	return [H2.getHamiltonians(mapping), geometry];


def main():
	mapping = raw_input('''Please specify the transformation ('BK' for Bravyi-Kitaev or 'JW' for Jordan-Wigner)\n''');

	while((mapping != 'BK') and (mapping != 'JW')):
		mapping = raw_input('''Could not understand input: {}. Please try again\n'''.format(mapping));

	atoms = int(raw_input('''Please specify the number of H atoms in your chain.\nNOTE: MUST BE AN INTEGER\n'''));

	while(atoms < 1):
		atoms = int(raw_input('''Could not understand input: {}. Please try again\n'''.format(atoms)));

	uniform = raw_input('''Please specify whether you would like the chain to be of uniform length (yes/no).\n''');

	while((uniform != "yes") and (uniform != "no")):
		uniform = raw_input('''Could not understand input: {}. Please try again\n'''.format(uniform));
	if(uniform == "yes"):
		uniform = True;
	else:
		uniform = False;

	bondLength = "random";
	if(uniform):
		bondLength = float(raw_input('''Please specify the bond length between atoms (angstroms).\nNOTE: MUST BE A NUMBER (DECIMAL VALUES ACCEPTED)\n'''));

		while(bondLength < 0):
			bondLength = float(raw_input('''Could not understand input: {}. Please try again\n'''.format(bondLength)));

	linear = raw_input('''Please specify whether you would like the chain to be straight (no y coordinates) (yes/no).\n''');

	while((linear != "yes") and (linear != "no")):
		linear = raw_input('''Could not understand input: {}. Please try again\n'''.format(linear));
	if(linear == "yes"):
		linear = True;
	else:
		linear = False;

	dimensions = raw_input('''Please specify whether you would like to model the molecule in 3 dimensions (yes) or 2 dimensions (no) (yes/no).\n''');

	while((dimensions != "yes") and (dimensions != "no")):
		dimensions = raw_input('''Could not understand input: {}. Please try again\n'''.format(dimensions));
	if(dimensions == "yes"):
		dimensions = True;
	else:
		dimensions = False;
	
	result = generateHChain(bondLength=bondLength, mapping=mapping, 
						   atoms=atoms, uniform=uniform, linear=linear, 
						   threedimensions=dimensions);

	print("\n\n ----- GEOMETRY ----- \n\n");
	print(result[1]);
	print("\n\n ----- HAMILTONIAN ----- \n\n");
	print(result[0]);

# main();


generateHChain(mapping="BK", atoms=1, linear=True, threedimensions=False);
















