#!/bin/bash

LMP=../../../../../src/lmp_serial
origdir=$(pwd)
for dn in nvst*/
do
    echo $dn
    cd $dn
    $LMP -i lammps.in > lammps.out
    cd $origdir
done
