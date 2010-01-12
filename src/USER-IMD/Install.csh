# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp -p fix_imd.cpp ..

  cp -p fix_imd.h ..

else if ($1 == 0) then

  rm ../fix_imd.cpp

  rm ../fix_imd.h

endif
