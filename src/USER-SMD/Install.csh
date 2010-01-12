# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp -p fix_smd.cpp ..

  cp -p fix_smd.h ..

else if ($1 == 0) then

  rm ../fix_smd.cpp

  rm ../fix_smd.h

endif
