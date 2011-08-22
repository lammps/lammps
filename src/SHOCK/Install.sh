# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_msst.cpp fix_nphug.cpp ..

  cp fix_msst.h fix_nphug.h ..

elif (test $1 = 0) then

  rm -f ../fix_msst.cpp ../fix_nphug.cpp

  rm -f ../fix_msst.h ../fix_nphug.h

fi
