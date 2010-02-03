# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p fix_smd.cpp ..

  cp -p fix_smd.h ..

elif (test $1 = 0) then

  rm ../fix_smd.cpp

  rm ../fix_smd.h

fi
