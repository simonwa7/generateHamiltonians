from generateCircuit import GenerateCircuit
import random

def createGeometry(bondLength, atoms, uniform, linear):
	index = 0;
	geometry = list();	

	#initialize the bondLength 
	if(bondLength == "random"):
		bondLength = random.uniform(0.1, 2);

	positionX = bondLength;
	positionY = 0;
	while(index < atoms):
		atom = tuple(("H", (positionX, positionY,0)));
		geometry.append(atom);
		index += 1;

		# Set the position for the next atom 
		if(uniform):
			positionX += bondLength;
		else:
			positionX += random.uniform(0.1, 2);
	print(geometry);


def generateHChain(bondLength=.7414, mapping="BK", atoms=2, uniform=True, linear=True):

	H2 = GenerateCircuit();
	H2.set_name("HChain");
	description = str(bondLength) + str(atoms) + str(uniform) + str(linear);
	H2.set_description(description);

	geometry = createGeometry(bondLength, atoms, uniform, linear);

	# print(bondLength);
	# H2.set_geometry([("H", (0, 0 ,0)), ("H", (bondLength, 0, 0))]);
	# H2.set_basis("sto-3g");
	# H2.set_multiplicity(1);
	# H2.set_charge(0);
	# H2.load_molecule();
	# H2.create_hamiltonians();
	# H2.create_circuits(mapping);
	# print(H2.getHamiltonians(mapping));

generateHChain(bondLength="random", atoms=5, uniform=False);