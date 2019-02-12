# LAMMPS LATTICE DYNAMICS COMMANDS

## THIRD ORDER TENSOR CALCULATOR

This directory contains the ingredients to calculate a third order tensor.  

Example:
```
$THIRD_ORDER=third_order #tensor output file
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in in.silicon -out out.silicon
combine.sh third_order
```

To test out a different silicon example:
```
$THIRD_ORDER=third_order
$LMP_FILE=amorphous_silicon.lmp
cp lmp_bank/$LMP_FILE ./silicon_input_file.lmp
NP=4 #number of processors
mpirun -np $NP lmp_mpi -in in.silicon -out out.silicon
bash combine.sh $THIRD_ORDER
```

## Requires: MANYBODY and MOLECULE packages
