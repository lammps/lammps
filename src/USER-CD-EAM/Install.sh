# Install/Uninstall package files in LAMMPS

if (test $1 = 1) then

  cp pair_cdeam.cpp ..

  cp pair_cdeam.h ..
  
elif (test $1 = 0) then

  rm -f ../pair_cdeam.cpp

  rm -f ../pair_cdeam.h

fi
