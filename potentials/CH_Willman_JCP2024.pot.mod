# DATE: 2024-08-09 UNITS: metal CONTRIBUTOR: Jonathan Tyler willman <willman@lanl.gov> CITATION: J.T. Willman, R. Perriot, C. Ticknor, "Atomic cluster expansion potential for large scale simulations of hydrocarbons under shock compression", The Journal of Chemical Physics, 161, 6, (2024). https://doi.org/10.1063/5.0213560

#STANDARD PACE
pair_style pace product chunksize 4096	
pair_coeff              * *  ./CH_ace_potential.yace C H 

