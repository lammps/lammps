# Install/unInstall package files in LAMMPS

if (test $1 == 1) then

  cp fix_srd.cpp ..
  cp fix_srd2.cpp ..
  cp fix_wall_srd.cpp ..

  cp fix_srd.h ..
  cp fix_srd2.h ..
  cp fix_wall_srd.h ..

elif (test $1 == 0) then

  rm ../fix_srd.cpp
  rm ../fix_srd2.cpp
  rm ../fix_wall_srd.cpp

  rm ../fix_srd.h
  rm ../fix_srd2.h
  rm ../fix_wall_srd.h

fi
