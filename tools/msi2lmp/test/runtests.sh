#!/bin/sh

BIOSYM_LIBRARY=../biosym_frc_files
MSI2LMP=../src/msi2lmp.exe
LAMMPS=../../../src/lmp_serial
CHECKDATA=./data-compare.pl

counter=0

# Class1 tests
${MSI2LMP} hydrogen-class1 -c 1 -p 2 	\
	|| counter=$(expr $counter + 1)
${LAMMPS} -log none -in in.hydrogen-class1	\
	|| counter=$(expr $counter + 1)
${CHECKDATA} hydrogen-class1.data reference/hydrogen-class1.data	\
	|| counter=$(expr $counter + 1)
${CHECKDATA} hydrogen-class1.data2 reference/hydrogen-class1.data2	\
	|| counter=$(expr $counter + 1)

# Class2 tests
${MSI2LMP} hydrogen-class2 -c 2 -p 2 -f compass_published	\
	|| counter=$(expr $counter + 1)
${LAMMPS} -log none -in in.hydrogen-class2	\
	|| counter=$(expr $counter + 1)
${CHECKDATA} hydrogen-class2.data reference/hydrogen-class2.data	\
	|| counter=$(expr $counter + 1)
${CHECKDATA} hydrogen-class2.data2 reference/hydrogen-class2.data2	\
	|| counter=$(expr $counter + 1)

echo "Total error count: $counter"
echo

# cleaning up.
rm log.* *.data *.data2
