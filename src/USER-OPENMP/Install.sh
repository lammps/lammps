# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p pair_omp.cpp ..

  cp -p pair_omp.h ..

elif (test $1 = 0) then

  rm ../pair_omp.cpp

  rm ../pair_omp.h

fi
