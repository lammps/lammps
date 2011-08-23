# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_msst.cpp ..
  cp fix_nphug.cpp ..

  cp fix_msst.h ..
  cp fix_nphug.h ..

elif (test $1 = 0) then

  rm -f ../fix_msst.cpp
  rm -f ../fix_nphug.cpp

  rm -f ../fix_msst.h
  rm -f ../fix_nphug.h

fi
