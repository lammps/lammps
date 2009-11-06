# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_dsmc.h ..

  cp pair_dsmc.cpp ..

  cp pair_dsmc.h ..

else if ($1 == 0) then

  rm ../style_dsmc.h
  touch ../style_dsmc.h

  rm ../pair_dsmc.cpp

  rm ../pair_dsmc.h

endif
