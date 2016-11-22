#!/bin/sh

MSI2LMP_LIBRARY=../frc_files
VALGRIND='valgrind -v --track-origins=yes --show-reachable=yes --leak-check=full'
MSI2LMP=../src/msi2lmp.exe
LAMMPS=../../../src/lmp_serial
CHECKDATA=./data-compare.pl

if [ ! -x $MSI2LMP ]
then
   echo "No executable $MSI2LMP"
   exit 1
fi

if [ ! -d $MSI2LMP_LIBRARY ]
then
   echo "No directory $MSI2LMP_LIBRARY"
   exit 1
fi

if [ ! -x $LAMMPS ]
then
   echo "No executable $LAMMPS"
   exit 1
fi

verbose=1
counter=0
errors=0

# Class1 tests with cvff
for m in hydrogen water h2-h2o ethane benzene naphthalene cnt-hexagonal crambin nylon phen3_cff97 hap_crystal
do \
    before=$errors
    vglog=${m}-class1.chk
    ${VALGRIND} --log-file=${vglog}				\
        ${MSI2LMP} ${m}-class1 -c 1 -p ${verbose}		\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class1		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class1.data reference/${m}-class1.data	\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class1.data2 reference/${m}-class1.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class1.data ${m}-class1.data2 log.${m}-class1
    leak=$(awk '/in use at exit:/ {num=$6;} END {print num;}' $vglog)
    [ $leak != 0 ] && echo "Memory still used: $leak"		\
        && errors=$(expr $errors + 1)
    viol=$(awk '/ERROR SUMMARY/ {num=$4;} END {print num;}' $vglog)
    [ $viol != 0 ] && echo "Valgrind errors: $viol"		\
        && errors=$(expr $errors + 1)
    [ $leak = 0 ] && [ $viol = 0 ] && rm ${vglog}
    counter=$(expr $counter + 6)
done

# Class1 tests with clayff
for m in PyAC_bulk
do \
    before=$errors
    vglog=${m}-clayff.chk
    ${VALGRIND} --log-file=${vglog}				\
        ${MSI2LMP} ${m}-clayff -c 1 -p ${verbose} -f clayff -n	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-clayff		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-clayff.data reference/${m}-clayff.data	\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-clayff.data2 reference/${m}-clayff.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-clayff.data ${m}-clayff.data2 log.${m}-clayff
    leak=$(awk '/in use at exit:/ {num=$6;} END {print num;}' $vglog)
    [ $leak != 0 ] && echo "Memory still used: $leak"		\
        && errors=$(expr $errors + 1)
    viol=$(awk '/ERROR SUMMARY/ {num=$4;} END {print num;}' $vglog)
    [ $viol != 0 ] && echo "Valgrind errors: $viol"		\
        && errors=$(expr $errors + 1)
    [ $leak = 0 ] && [ $viol = 0 ] && rm ${vglog}
    counter=$(expr $counter + 6)
done

# OPLS-AA tests 
for m in ethane decane
do \
    before=$errors
    vglog=${m}-oplsaa.chk
    ${VALGRIND} --log-file=${vglog}				\
        ${MSI2LMP} ${m}-oplsaa -c 0 -p ${verbose} -f oplsaa -n	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-oplsaa		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-oplsaa.data reference/${m}-oplsaa.data	\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-oplsaa.data2 reference/${m}-oplsaa.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-oplsaa.data ${m}-oplsaa.data2 log.${m}-oplsaa
    leak=$(awk '/in use at exit:/ {num=$6;} END {print num;}' $vglog)
    [ $leak != 0 ] && echo "Memory still used: $leak"		\
        && errors=$(expr $errors + 1)
    viol=$(awk '/ERROR SUMMARY/ {num=$4;} END {print num;}' $vglog)
    [ $viol != 0 ] && echo "Valgrind errors: $viol"		\
        && errors=$(expr $errors + 1)
    [ $leak = 0 ] && [ $viol = 0 ] && rm ${vglog}
    counter=$(expr $counter + 6)
done

# Class2 tests with compass
for m in hydrogen ethane benzene naphthalene cnt-hexagonal
do \
    before=$errors
    vglog=${m}-class2a.chk
    ${VALGRIND} --log-file=${vglog}				\
        ${MSI2LMP} ${m}-class2a -c 2 -p ${verbose} -f compass_published	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class2a		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2a.data reference/${m}-class2a.data		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2a.data2 reference/${m}-class2a.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class2a.data ${m}-class2a.data2 log.${m}-class2a
    leak=$(awk '/in use at exit:/ {num=$6;} END {print num;}' $vglog)
    [ $leak != 0 ] && echo "Memory still used: $leak"		\
        && errors=$(expr $errors + 1)
    viol=$(awk '/ERROR SUMMARY/ {num=$4;} END {print num;}' $vglog)
    [ $viol != 0 ] && echo "Valgrind errors: $viol"		\
        && errors=$(expr $errors + 1)
    [ $leak = 0 ] && [ $viol = 0 ] && rm ${vglog}
    counter=$(expr $counter + 6)
done

# Class2 tests with pcff
for m in water h2-h2o ethane benzene naphthalene cnt-hexagonal hap_crystal
do \
    before=$errors
    vglog=${m}-class2b.chk
    ${VALGRIND} --log-file=${vglog}				\
        ${MSI2LMP} ${m}-class2b -c 2 -p ${verbose} -f pcff	\
        || errors=$(expr $errors + 1)
    ${LAMMPS} -log none -screen none -in in.${m}-class2b		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2b.data reference/${m}-class2b.data		\
        || errors=$(expr $errors + 1)
    ${CHECKDATA} ${m}-class2b.data2 reference/${m}-class2b.data2	\
        || errors=$(expr $errors + 1)
    [ $before -eq $errors ] && rm ${m}-class2b.data ${m}-class2b.data2 log.${m}-class2b
    leak=$(awk '/in use at exit:/ {num=$6;} END {print num;}' $vglog)
    [ $leak != 0 ] && echo "Memory still used: $leak"		\
        && errors=$(expr $errors + 1)
    viol=$(awk '/ERROR SUMMARY/ {num=$4;} END {print num;}' $vglog)
    [ $viol != 0 ] && echo "Valgrind errors: $viol"		\
        && errors=$(expr $errors + 1)
    [ $leak = 0 ] && [ $viol = 0 ] && rm ${vglog}
    counter=$(expr $counter + 6)
done

echo "Total error count: $errors / $counter"
echo

echo "Total error count: $errors / $counter"
echo

