#!/bin/sh

BIOSYM_LIBRARY=../biosym_frc_files
MSI2LMP=../src/msi2lmp.exe
LAMMPS=../../../src/lmp_serial
CHECKDATA=./data-compare.pl

counter=0

# Class1 tests
for m in hydrogen water
do \
    before=$counter
    ${MSI2LMP} ${m}-class1 -c 1 -p 2 	\
        || counter=$(expr $counter + 1)
    ${LAMMPS} -log none -in in.${m}-class1	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class1.data reference/${m}-class1.data	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class1.data2 reference/${m}-class1.data2	\
        || counter=$(expr $counter + 1)
    [ $before -eq $counter ] && rm ${m}-class1.data ${m}-class1.data2 log.${m}-class1
done

# Class2 tests with compass
for m in hydrogen
do \
    before=$counter
    ${MSI2LMP} ${m}-class2 -c 2 -p 2 -f compass_published	\
        || counter=$(expr $counter + 1)
    ${LAMMPS} -log none -in in.${m}-class2	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class2.data reference/${m}-class2.data	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class2.data2 reference/${m}-class2.data2	\
        || counter=$(expr $counter + 1)
    [ $before -eq $counter ] && rm ${m}-class2.data ${m}-class2.data2 log.${m}-class2
done

# Class2 tests with cff91
for m in water
do \
    before=$counter
    ${MSI2LMP} ${m}-class2 -c 2 -p 2 -f cff91	\
        || counter=$(expr $counter + 1)
    ${LAMMPS} -log none -in in.${m}-class2	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class2.data reference/${m}-class2.data	\
        || counter=$(expr $counter + 1)
    ${CHECKDATA} ${m}-class2.data2 reference/${m}-class2.data2	\
        || counter=$(expr $counter + 1)
    [ $before -eq $counter ] && rm ${m}-class2.data ${m}-class2.data2 log.${m}-class2
done

echo "Total error count: $counter"
echo

