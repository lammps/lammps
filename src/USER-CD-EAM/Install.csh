# Install/Uninstall package classes in LAMMPS

if ($1 == 1) then

  cp style_user_cd_eam.h ..
  
  cp pair_cdeam.cpp ..
  cp pair_cdeam.h ..
  
else if ($1 == 0) then

  rm ../style_user_cd_eam.h
  touch ../style_user_cd_eam.h

  rm ../pair_cdeam.cpp
  rm ../pair_cdeam.h

endif
