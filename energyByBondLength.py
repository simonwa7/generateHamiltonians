from generateCircuit import GenerateCircuit
import random
import math
from openfermion.transforms import get_sparse_operator
from openfermion.utils import get_ground_state

MIN_LENGTH = 0.2;
# MAX_LENGTH = 1.4828;
MAX_LENGTH = 4;
ITERATIONS = 10;

MIN_ANGLE = 0;
MAX_ANGLE = 10; # (theta-1)*10 = degrees

BASIS = '3-21g'

def getLowest(geometry, middleLength):
	molecule = GenerateCircuit();
	molecule.set_name("H2+H");
	molecule.set_geometry(geometry);
	molecule.set_basis(BASIS);
	molecule.set_multiplicity(2);
	molecule.set_charge(0);
	molecule.load_molecule();

	energy = molecule.molecule.hf_energy;
	molecule.create_hamiltonians();

	sparse = get_sparse_operator(molecule.fermion_hamiltonian)
	values = get_ground_state(sparse)[0]

	return [values, energy];

def linearH3():
	results = [];

	middleLength = MIN_LENGTH;
	totalLength = MAX_LENGTH;

	while(middleLength < totalLength):
		geometry = createGeometryLinear(totalLength, middleLength);
		data = getLowest(geometry, middleLength);
		result = [middleLength, data[0], data[1]]
		results.append(result);
		middleLength += (MAX_LENGTH-MIN_LENGTH)/(ITERATIONS);

	file = open("linearEnergies.txt", "w");
	for data in results:
		file.write("{} {} {}\n".format(data[0], data[1], data[2]));
	file.close();

def rotatingH3():
	results = [];

	bondLength = 0.7414; # middle atom at normal bond length for H2

	for theta in range(MIN_ANGLE, MAX_ANGLE):
		theta *= 10;
		radian = theta*math.pi/180;
		geometry = createGeometryRotated(bondLength, radian);
		data = getLowest(geometry, bondLength, radian=radian);
		result = [theta, data[0], data[1]];
		results.append(result);

	file = open("rotatedEnergies.txt", "w");
	for data in results: 
		file.write("{} {} {}\n".format(data[0], data[1], data[2]));
	file.close();

def wobblingH3():
	middleLength = MIN_LENGTH;
	totalLength = MAX_LENGTH;

	dataFilename = "GroundStateEnergies_" + BASIS + "_" + MIN_LENGTH + "_" + MAX_LENGTH + "_" + MIN_ANGLE + "_" + MAX_ANGLE + "_" + ITERATIONS;
	file = open(dataFilename, "w");
	
	while(middleLength < totalLength-MIN_LENGTH):
		bondLength = totalLength-middleLength;
		
		for theta in range(MIN_ANGLE, MAX_ANGLE):
			print(middleLength, theta);

			geometry = list();
			geometry.append(("H", (0.,0.,0.)));
			geometry.append(("H", (0.,0.,middleLength)));
			
			radian = theta*10*math.pi/180;
			positionX = middleLength + (bondLength*math.cos(radian));
			positionY = bondLength*math.sin(radian);
			geometry.append(("H", (0., positionY, positionX)));
			print(geometry)

			hf_energy, ground_state  = getLowest(geometry, middleLength);

			result = [middleLength, theta, data[0], data[1]];
			file.write("{} {} {} {}\n".format(middleLength, str(theta*10), hf_energy, ground_state));

		middleLength += (MAX_LENGTH-(2*MIN_LENGTH))/(ITERATIONS)
	file.close();

# linearH3();
# rotatingH3();
wobblingH3();