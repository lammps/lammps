mpirun -np 4 $LMP -in in.langevin.metal -p 4x1 -log log.langevin.metal -screen screen
mpirun -np 4 $LMP -in in.pimd-langevin.metal -p 4x1 -log log.pimd-langevin.metal -screen screen

