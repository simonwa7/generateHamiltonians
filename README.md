Author: William A. Simon

Date: 11/7/2018
# Generate Hamiltonians
### Purpose: 
The purpose of this script is to easily generate qubit hamiltonians for varying molecules made solely of Hydrogen chains. 

#### Usage:
```
python generateHamiltonians.py
```

Through the command line, a client can define various parameters such as:

* transformation - Bravyi-Kitaev or Jordan-Wigner</li>
* number of hydrogen atoms - integer</li>
* whether or not the molecules are separated by a uniform distance - yes/no</li>
* if they are separated by a uniform distance, what is that bond length - angstroms (float)</li>
* if the chain is straight or randomly bent - yes/no</li>

#### Examples:
* [2 Dimensions or Linear Chains](https://github.com/simonwa7/generateHamiltonians/blob/master/exampleOutput.txt "2D Modeling")
* [3 Dimensions](https://github.com/simonwa7/generateHamiltonians/blob/master/exampleOutput3D.txt "3D Modeling")

#### Helpful Tools:
* Check out [this online program](https://www.geogebra.org/3d?lang=en "3D Graphing") to take a look at what your molecule looks like!

#### Dependencies:
This script utilizes the open source projects OpenFerimon and OpenFermion-psi4 to generate molecular integral data and perform the necessary transformations.

#### Installation: 
with python and pip installed, simply run:
```
pip install -r requirements.txt
```

#### TODO:
