from generateCircuit import GenerateCircuit
import random
import math

# minimum bond angle is pi-MIN_THETA
# maximum bond angle is always 180
MIN_THETA = -math.pi/3;
MAX_THETA = math.pi/3;

def generateHChain(bondLength=.7414, mapping="BK", atoms=2, uniform=True, 
				   linear=True, theta=0):

	H2 = GenerateCircuit();
	H2.set_name("HChain");
	description = str(bondLength) + str(atoms) + str(uniform) + str(linear);
	H2.set_description(description);

	geometry = [("H", (0,0,0)), ("H", (bondLength, 0, 0)), ("H", (bondLength+(bondLength*math.cos(theta)), bondLength*math.sin(theta), 0))];
	file = open("geometries.txt", "a");
	file.write("{} {}\n".format(theta, geometry[2]));
	file.close();
	H2.set_geometry(geometry);
	H2.set_basis("sto-3g");
	H2.set_multiplicity(1);
	H2.set_charge(atoms-2);
	H2.load_molecule();
	H2.create_hamiltonians();
	H2.create_circuits(mapping);
	# print(H2.getHamiltonians(mapping));
	# print("\n\n");
	values, vectors = H2.getEigen(mapping);
	print("\n\nEIGENVALUES FOR ANGLE {} Degrees".format(theta*180/math.pi))
	print(values)
	# print("\nEIGENVECTORS")
	# print(vectors)
	return None;


def main():
	result = [];
	angles = [];
	for theta in range(0, 19):
		theta *= 10;
		radian = theta*math.pi/180;
		# angles.append(radian);
		result.append(generateHChain(mapping="JW", atoms=3, linear=False, theta=radian));

	# print("\n\n ----- EIGENVALUES ----- \n\n");
	# print(result);

	# file = open("data.txt", "a");
	# for index, angle in enumerate(angles):
	# 	file.write("{} {}\n".format(angles[index], result[index][0]));
	# file.close();
main();
