from generateCircuit import GenerateCircuit
import random
import math

MIN_LENGTH = 0.1;
MAX_LENGTH = 2;

MIN_THETA = 0;
MAX_THETA = 2*math.pi;

def createGeometry(bondLength, atoms, uniform, linear):
	index = 1;
	geometry = list();	

	# Add first H atom at 0, 0, 0
	geometry.append(tuple(("H", (0,0,0))));

	#initialize the bondLength 
	if(bondLength == "random"):
		bondLength = random.uniform(MIN_LENGTH, MAX_LENGTH);

	positionX = 0;
	positionY = 0;
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
				positionX += bondLength*math.cos(theta);
				positionY += bondLength*math.sin(theta);
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
				positionX += bondLength*math.cos(theta);
				positionY += bondLength*math.sin(theta);

		# create atom tuple and push it to the list
		atom = tuple(("H", (positionX, positionY,0)));
		geometry.append(atom);
		index += 1;


	return geometry;


def generateHChain(bondLength=.7414, mapping="BK", atoms=2, uniform=True, linear=True):

	H2 = GenerateCircuit();
	H2.set_name("HChain");
	description = str(bondLength) + str(atoms) + str(uniform) + str(linear);
	H2.set_description(description);

	geometry = createGeometry(bondLength, atoms, uniform, linear);

	H2.set_geometry(geometry);
	H2.set_basis("sto-3g");
	H2.set_multiplicity(1);
	H2.set_charge(atoms-2);
	H2.load_molecule();
	H2.create_hamiltonians();
	H2.create_circuits(mapping);
	print(H2.getHamiltonians(mapping));

generateHChain(bondLength=1, mapping="BK", atoms=5, uniform=False, linear=False);




