# Install/Uninstall package classes in LAMMPS

if ($1 == 1) then

  cp pair_cdeam.cpp ..
  cp pair_cdeam.h ..
  
else if ($1 == 0) then

  rm ../pair_cdeam.cpp
  rm ../pair_cdeam.h

endif
