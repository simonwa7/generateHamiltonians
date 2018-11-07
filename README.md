Author: William A. Simon
Date: 11/7/2018
<h1>Generate Hamiltonians</h1>
<h3>
Purpose: The purpose of this script is to easily generate qubit hamiltonians for varying molecules made solely of Hydrogen chains. </h3>

<h4>Usage: Through the command line, a client can define various parameters such as:</h4>
<ul>
	<li>transformation - Bravyi-Kitaev or Jordan-Wigner</li>
	<li>number of hydrogen atoms - integer</li>
	<li>whether or not the molecules are separated by a uniform distance - yes/no</li>
	<li>if they are separated by a uniform distance, what is that bond length - angstroms (float)</li>
	<li>if the chain is straight or randomly bent - yes/no</li>
</ul>

<h4>Dependencies: This script utilizes the open source projects OpenFerimon and OpenFermion-psi4 to generate molecular integral data and perform the necessary transformations.</h4>

<h4>TODO:</h4>
