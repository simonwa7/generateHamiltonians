from generateCircuit import GenerateCircuit
import random
import math

MIN_LENGTH = 0.2;
# MAX_LENGTH = 1.4828;
MAX_LENGTH = 4;
ITERATIONS = 25;

MIN_ANGLE = 0;
MAX_ANGLE = 10; # (theta-1)*10 = degrees

def createGeometryLinear(totalLength, middleLength): 
	geometry = list();	

	# Add first H atom at 0, 0, 0
	geometry.append(tuple(("H", (0.,0.,0.))));
	geometry.append(tuple(("H", (0., 0., middleLength))));
	geometry.append(tuple(("H", (0., 0., totalLength))));

	return geometry;

def createGeometryRotated(bondLength, radian):
	geometry = list();

	geometry.append(tuple(("H", (0.,0.,0.))));
	geometry.append(tuple(("H", (0., 0., bondLength))));
	geometry.append(tuple(("H", (0., bondLength*math.sin(radian), 
								 bondLength+(bondLength*math.cos(radian))))));
	return geometry;

def getLowest(geometry, middleLength, radian=0):
	molecule = GenerateCircuit();
	molecule.set_name("oscillatingHChain");
	description = "wobbly_middlelocation="+str(middleLength)+"charge="+str(-1)+"angle="+str(radian);
	molecule.set_description(description);
	molecule.set_geometry(geometry);
	molecule.set_basis("sto-3g");
	molecule.set_multiplicity(1);
	molecule.set_charge(-1);
	molecule.load_molecule();
	energy = molecule.molecule.hf_energy;
	molecule.create_hamiltonians();
	molecule.create_circuits("BK");
	values = molecule.getLowestEigen("BK");
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
	results = [];

	middleLength = MIN_LENGTH;
	totalLength = MAX_LENGTH;

	
	while(middleLength < totalLength):
		bondLength = totalLength-middleLength;
		for theta in range(MIN_ANGLE, MAX_ANGLE):
			geometry = list();
			geometry.append(("H", (0.,0.,0.)));
			geometry.append(("H", (0.,0.,middleLength)));
			print(middleLength, theta);
			theta *= 10;
			radian = theta*math.pi/180;
			positionX = middleLength + (bondLength*math.cos(radian));
			positionY = bondLength*math.sin(radian);
			geometry.append(("H", (0., positionY, positionX)));
			print(geometry)
			data = getLowest(geometry, middleLength, radian=radian);
			# result = [0, middleLength, positionY, positionX, theta, data[0], data[1]];
			result = [middleLength, theta, data[0], data[1]];
			results.append(result);

		middleLength += (MAX_LENGTH-MIN_LENGTH)/(ITERATIONS)

	file = open("wobblingEnergies.txt", "w");
	for data in results: 
		file.write("{} {} {} {}\n".format(data[0], data[1], data[2], data[3]));
	file.close();

# linearH3();
# rotatingH3();
wobblingH3();