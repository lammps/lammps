# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_xtc.h ..

  cp dump_xtc.cpp ..

  cp dump_xtc.h ..

else if ($1 == 0) then

  rm ../style_xtc.h
  touch ../style_xtc.h

  rm ../dump_xtc.cpp

  rm ../dump_xtc.h

endif
