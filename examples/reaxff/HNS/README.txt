# Pure HNS crystal, ReaxFF tests for benchmarking LAMMPS

Authors: Mitchell Wood, Stan Moore, and Aidan Thompson, Sandia National Labs
Date: 2017-10-17
Questions: Mitchell Wood, mitwood@sandia.gov

---------ReaxFF Benchmark Files---------
1) in.reaxc.hns
2) data.hns-equil
3) ffield.reax.hns

---------File Description---------

1) LAMMPS input file for this test. 
  The file is read line by line looking for keywords to set up this run. 
  It will read in the configuration given by the argument of the read_data command, which is supplied in this distribution.
  The type of simulation is set by the 'fix' commands, dynamic charges are controlled with 'fix qeq' and the integration style is given as 'fix nve' here.
  More information about each of the individual commands can be found online at www.lammps.org in the user manual section.

  *There are four free variables in this file, three of which control the size of the simulation and the last will dictate how many MD time steps are taken.
  *The size of the system is controlled by the 'replicate' command given the values of $x, $y and $z.
  *The number of timesteps taken is controlled by the 'run' command given the value of $t

  It is worth noting that these four free variables can be set at the command line when the simulation is invoked rather than editing the file by hand prior to each run.

  Example syntax:
    lmp_serial -in in.reaxc.hns -v x 2 -v y 2 -v z 2 -v t 100

2) LAMMPS Data file for crystalline HNS
  This file matches the LAMMPS data format, more information about this data structure can be found at www.lammps.org
  
  This particular data file is of the energetic material Hexanitrostilbene (HNS) with atom_style charge (id type q x y z).
  The file contains eight molecules (2 unit cells).

    Chemical Name: Hexanitrostilbene
    Molecule Composition:  C14H6N6O12
    IUPAC Name: 1,3,5-Trinitro-2-[2-(2,4,6-trinitrophenyl)ethenyl]benzene    
    
    Phyical Properties (Gerard F., Hardy A. Acta Cryst. (1988) 1283-1287)
    Density: 1.741 g/cc
    Crystal Structure: Monoclinic P2(1)/c
    Molecules per Unit Cell: 4
    Lattice Constants: a=22.326
           b=5.5706
           c=14.667
           beta=110.04 deg

3) ReaxFF force field file.
  Details about this particular parameterization can be found in T.-R. Shan, R. R. Wixom, A. P. Thompson, "Atomistic Simulation of Nanoscale Void-Enhanced Initiation in Hexanitrostilbene", Proc. 15th International Detonation Symposium, pp. 962-969, SAND Number: SAND2014-15518J
