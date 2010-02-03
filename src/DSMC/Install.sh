# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp pair_dsmc.cpp ..

  cp pair_dsmc.h ..

elif (test $1 = 0) then

  rm ../pair_dsmc.cpp

  rm ../pair_dsmc.h

fi
