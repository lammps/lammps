# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_srd.cpp ..
  cp fix_wall_srd.cpp ..

  cp fix_srd.h ..
  cp fix_wall_srd.h ..

elif (test $1 = 0) then

  rm -f ../fix_srd.cpp
  rm -f ../fix_wall_srd.cpp

  rm -f ../fix_srd.h
  rm -f ../fix_wall_srd.h

fi
