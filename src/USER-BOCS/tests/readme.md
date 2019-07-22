This folder holds input and output file sets useful for testing user-bocs code.

In the absence of automated unit testing, having test benchmarks is the next best thing.

Folder / description:

###cubic_spline-file_test
The base file set was supplied by Michael Delyser as a demo input file for testing
cubic spline basis analysis.
- cgheptane.lmp
- data.cgheptane
- delta_Fv.dat
- lammps_ang_CT-CM-CT.table
- lammps_bond_CT-CM.table
- lammps_CM_CM.table
- lammps_CT_CM.table
- lammps_CM_CT.table
- cgheptane.baseFilesetOutput.txt -- sample output captured 20190722 from master branch code

In addition, are now a set of either stand-alone .lmp files or pairs of .lmp/.dat files that
exercise various bits of error handling code. I'm hoping filenames are self-explanatory. Look
for .bad in the filename.

