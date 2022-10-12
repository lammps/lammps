#!/bin/bash

## run REMD using LAMMPS
mpirun -np 16 ~/mysoftware/lammps/src/lmp_mpi -partition 16x1 -in in.peptide -log log.peptide

## collect all energies from different replica logs
echo ; echo
echo "Parsing energies from replica logs"
python parse_ene.py temps.txt log.peptide

## run the reordering tool to get reordered trajectories @ 200 K, 276 K, 400 K
echo ; echo
mpirun -np 16 python ../reorder_remd_traj.py peptide -logfn log.peptide -tfn temps.txt -ns 10 -nw 20 -np 1000 -ot 200 276 400 -logw -e ene.peptide -od ./output
