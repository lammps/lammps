# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  cp fix_phonon.h ..
  cp fix_phonon.cpp ..

elif (test $1 = 0) then

  rm -f ../fix_phonon.h
  rm -f ../fix_phonon.cpp

fi
