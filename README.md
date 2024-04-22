# UNIPAR_V5
UNIPAR_V5 submodule with SOA_H5 and SOA_AE3

This directory contains the scource code for UNIPAR_V5 including:
	1. UNIPAR_V5 (Main module which will call other submodules)
	2. SOA_H5 (Submodule to simulate aerosol phase reaction, polymerization in organic phase and acid-catalyzed polymerization in aqueous phase)
	3. SOA_AE3 	(Submodule to simulate gas-organic partitioning of organic species using Newtonian model)

Main module (UNIPAR_V5) will call paramters for a specific precursor found in the SOA_PAR folder (i.e. SOA_PAR/SOAPAR_AU15)
	The SOA_PAR file contains paramters for each lumping species:
	1. Molecular weight
	2. O:C ratio
	3. Hydrogen Bonding
	4. NOx and aging dependent stoicometric coefficient equation
	
Main module will also call for an input file for a specific chamber simulation data found in the input folder (i.e. input/MMDDYY_East)
	which contains concentrations of HO2 (PPB), RO2 (PPB), change in ROG (dROG, PPB), Temp (K), RH (0-1), dSO4 (umol m^(-3)), dNH4 (umol m^(-3))
	Each column has a value for each 6 minute cycle. The preset is 100 cycles. 
