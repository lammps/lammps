# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp pair_dsmc.cpp ..

  cp pair_dsmc.h ..

else if ($1 == 0) then

  rm ../pair_dsmc.cpp

  rm ../pair_dsmc.h

endif
