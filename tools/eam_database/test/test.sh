#!/usr/bin/env bash

elements=(Cu Ag Au Ni Pd Pt Al Pb Fe Mo Ta W Mg Co Ti Zr Cr)
LMP=${LMP-../../../src/lmp_serial}

rm -r tmp
for ((i=0; i< ${#elements[@]}; i=$i+1))
do
    for (( j=$i; j<${#elements[@]}; j=$j+1))
    do
        e1=${elements[${i}]}
        e2=${elements[${j}]}
        echo "e1 = ${e1} e2 = ${e2}"
        pythonfile="${e1}${e2}.eam.alloy.python"
        fortranfile="${e1}${e2}.eam.alloy.fortran"
        mkdir tmp;
        cp eamDatabase.py create_eam.py in.lmp data.lmp a.out EAM.input EAM_code tmp/
        cd tmp
        sed -i "s/@ELEM1@/${e1}/g" EAM.input
        sed -i "s/@ELEM2@/${e2}/g" EAM.input
        ./a.out < EAM.input
        mv "${e1}${e2}.eam.alloy" ${fortranfile}
        python create_eam.py -n ${e1} ${e2}
        mv "${e1}${e2}.eam.alloy" ${pythonfile}
        sed -i "s/@ELEM1@/${e1}/g" in.lmp
        sed -i "s/@ELEM2@/${e2}/g" in.lmp
        sed -i "s/@PYTHONFILE@/${pythonfile}/g" in.lmp
        sed -i "s/@FORTRANFILE@/${fortranfile}/g" in.lmp
        ${LMP} -i in.lmp
        cd ../
        rm -r tmp
    done
done

