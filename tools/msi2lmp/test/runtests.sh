#!/bin/sh

BIOSYM_LIBRARY=../biosym_frc_files
MSI2LMP=../src/msi2lmp.exe
LAMMPS=../../../src/lmp_serial
CHECKDATA=./data-compare.pl

verbose=1
counter=0
errors=0

# Class1 tests with cvff
for m in hydrogen water h2-h2o ethane benzene naphthalene crambin nylon phen3_cff97
do \
    before=$errors
    ${MSI2LMP} ${m}-class1 -c 1 -p ${verbose}			\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class1		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class1.data reference/${m}-class1.data	\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class1.data2 reference/${m}-class1.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class1.data ${m}-class1.data2 log.${m}-class1
    counter=$(expr $counter + 4)
done

# Class1 tests with clayff
for m in PyAC_bulk
do \
    before=$errors
    ${MSI2LMP} ${m}-clayff -c 1 -p ${verbose} -f clayff	-n	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-clayff		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-clayff.data reference/${m}-clayff.data	\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-clayff.data2 reference/${m}-clayff.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-clayff.data ${m}-clayff.data2 log.${m}-clayff
    counter=$(expr $counter + 4)
done

# Class2 tests with compass
for m in hydrogen ethane benzene naphthalene
do \
    before=$errors
    ${MSI2LMP} ${m}-class2a -c 2 -p ${verbose} -f compass_published	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class2a		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2a.data reference/${m}-class2a.data		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2a.data2 reference/${m}-class2a.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class2a.data ${m}-class2a.data2 log.${m}-class2a
    counter=$(expr $counter + 4)
done

# Class2 tests with pcff
for m in water h2-h2o ethane benzene naphthalene
do \
    before=$errors
    ${MSI2LMP} ${m}-class2b -c 2 -p ${verbose} -f pcff			\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class2b		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2b.data reference/${m}-class2b.data		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2b.data2 reference/${m}-class2b.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class2b.data ${m}-class2b.data2 log.${m}-class2b
    counter=$(expr $counter + 4)
done

echo "Total error count: $errors / $counter"
echo

