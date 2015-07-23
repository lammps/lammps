#This script implements the BKS pair potential for various silicon dioxide compounds. Inner part is fixed with a harmonic potential. Long range Coulomb interactions are evaluated with the pppm method.

#Pair Potentials
pair_style		hybrid/overlay buck/coul/long ${cut_off} table linear 39901
pair_coeff		1 1 buck/coul/long 0.0 1.0 0.0								#No interactions between Si atoms 
pair_coeff		1 2 buck/coul/long 18003.757200 0.205205 133.538100	
pair_coeff		2 2 buck/coul/long 1388.773000  0.362319 175.000000					#BKS interaction in PRL 64 1955 (1990)
pair_modify		shift yes
pair_coeff		1 2 table potential_SiO2.TPF Si-O ${cut_off}
pair_coeff		2 2 table potential_SiO2.TPF O-O ${cut_off}						#See the potential file for more information
kspace_style		pppm 1.0e-4 

#Neighbor style
neighbor		2.0 bin
neigh_modify		check yes every 1 delay 0 page 100000 one 2000
