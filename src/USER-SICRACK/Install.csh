# Install/unInstall package "USER-SICRACK" classes in LAMMPS

if (test $1 == 1) then

  cp fix_sicrack.cpp ..
  cp fix_sicrack.h ..

elif (test $1 == 0) then

  rm ../fix_sicrack.cpp
  rm ../fix_sicrack.h

fi
