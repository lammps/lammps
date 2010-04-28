# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p pair_omp.cpp ..
  cp -p pair_lj_cut_omp.cpp ..

  cp -p pair_omp.h ..
  cp -p pair_lj_cut_omp.h ..

elif (test $1 = 0) then

  rm ../pair_omp.cpp
  rm ../pair_lj_cut_omp.cpp

  rm ../pair_omp.h
  rm ../pair_lj_cut_omp.h

fi
