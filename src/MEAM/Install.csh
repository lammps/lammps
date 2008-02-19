# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_meam.h ..

  cp pair_meam.cpp ..

  cp pair_meam.h ..

else if ($1 == 0) then

  rm ../style_meam.h
  touch ../style_meam.h

  rm ../pair_meam.cpp

  rm ../pair_meam.h

endif
