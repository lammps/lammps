# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp -p style_user_smd.h ..

  cp -p fix_smd.cpp ..

  cp -p fix_smd.h ..

else if ($1 == 0) then

  rm ../style_user_smd.h
  touch ../style_user_smd.h

  rm ../fix_smd.cpp

  rm ../fix_smd.h

endif
