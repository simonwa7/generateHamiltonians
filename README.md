Author: William A. Simon<br>
Date: 11/7/2018
<h1>Generate Hamiltonians</h1>
<h3>Purpose:</h3> The purpose of this script is to easily generate qubit hamiltonians for varying molecules made solely of Hydrogen chains. 

<h4>Usage:</h4> 
<br>
`python generateHamiltonians.py`
<br>
Through the command line, a client can define various parameters such as:

<h6><ul>
	<li>transformation - Bravyi-Kitaev or Jordan-Wigner</li>
	<li>number of hydrogen atoms - integer</li>
	<li>whether or not the molecules are separated by a uniform distance - yes/no</li>
	<li>if they are separated by a uniform distance, what is that bond length - angstroms (float)</li>
	<li>if the chain is straight or randomly bent - yes/no</li>
</ul></h6>

<h4>Dependencies:</h4> This script utilizes the open source projects OpenFerimon and OpenFermion-psi4 to generate molecular integral data and perform the necessary transformations.

<h4>Installation: </h4>with python and pip installed, simply run:
`pip install -r requirements.txt`

<h4>TODO:</h4>
