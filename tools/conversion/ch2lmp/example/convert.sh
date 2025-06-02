#!/bin/sh
rm -f 1ac7.restart* 1ac7.dump
../charmm2lammps.pl -water -ions -border=2.0 all27_na 1ac7
