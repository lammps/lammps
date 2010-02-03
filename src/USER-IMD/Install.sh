# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p fix_imd.cpp ..

  cp -p fix_imd.h ..

elif (test $1 = 0) then

  rm ../fix_imd.cpp

  rm ../fix_imd.h

fi
