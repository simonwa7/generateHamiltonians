from OpenFermionWrapper import OpenFermionWrapper
import math

MIN_LENGTH = 0.2
MAX_LENGTH = 3
ITERATIONS = 10

MIN_ANGLE = -16 # (theta-1)*10 = degrees
MAX_ANGLE = 18

BASIS = 'sto-3g'

def getEnergy(geometry):
	molecule = OpenFermionWrapper()
	molecule.set_name("H2+H")
	molecule.set_geometry(geometry)
	molecule.set_basis(BASIS)
	molecule.set_multiplicity(2)
	molecule.set_charge(0)
	molecule.load_molecule()

	energy = molecule.molecule.hf_energy
	molecule.create_hamiltonians()

	molecule.set_ground_state_energy()
	eigenvalue = molecule.set_ground_state_energy()
	return [energy, eigenvalue]

def main():
	middleLength = MIN_LENGTH
	totalLength = MAX_LENGTH

	dataFilename = "GroundStateEnergies_" + BASIS + "_" + str(MIN_LENGTH) + "_"
	dataFilename +=	str(MAX_LENGTH) + "_" + str(MIN_ANGLE) + "_"
	dataFilename += str(MAX_ANGLE) + "_" + str(ITERATIONS)

	file = open(dataFilename, "w")
	
	while(middleLength < totalLength-MIN_LENGTH):
		bondLength = totalLength-middleLength
		
		for theta in range(MIN_ANGLE, MAX_ANGLE):
			print(middleLength, theta)

			geometry = list()
			geometry.append(("H", (0.,0.,0.)))
			geometry.append(("H", (0.,0.,middleLength)))
			
			radian = theta*10*math.pi/180
			positionX = middleLength + (bondLength*math.cos(radian))
			positionY = bondLength*math.sin(radian)
			geometry.append(("H", (0., positionY, positionX)))
			print(geometry)

			hf_energy, ground_state  = getEnergy(geometry)

			file.write("{} {} {} {}\n".format(middleLength, str(theta*10), 
											  hf_energy,    ground_state));

		middleLength += (MAX_LENGTH-(2*MIN_LENGTH))/(ITERATIONS)
	file.close();

main()