LMP_EXEC=lmp_electrode
tags=(p3m ew3dc ew2d)

for t in ${tags[@]}
  do
    export TAG=$t
    mpirun -np 4 $LMP_EXEC -i conp.in
  done

paste printout_*
