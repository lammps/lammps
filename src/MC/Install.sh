# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp pair_dsmc.cpp ..

  cp pair_dsmc.h ..

elif (test $1 = 0) then

  rm -f ../pair_dsmc.cpp

  rm -f ../pair_dsmc.h

fi
