# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_msst.cpp ..

  cp fix_msst.h ..

elif (test $1 = 0) then

  rm ../fix_msst.cpp

  rm ../fix_msst.h

fi
