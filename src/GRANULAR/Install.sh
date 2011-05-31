# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_freeze.cpp ..
  cp fix_pour.cpp ..
  cp fix_wall_gran.cpp ..
  cp pair_gran_hertz_history.cpp ..
  cp pair_gran_hooke.cpp ..
  cp pair_gran_hooke_history.cpp ..

  cp fix_freeze.h ..
  cp fix_pour.h ..
  cp fix_wall_gran.h ..
  cp pair_gran_hertz_history.h ..
  cp pair_gran_hooke.h ..
  cp pair_gran_hooke_history.h ..

elif (test $1 = 0) then

  rm -f ../fix_freeze.cpp
  rm -f ../fix_pour.cpp
  rm -f ../fix_wall_gran.cpp
  rm -f ../pair_gran_hertz_history.cpp
  rm -f ../pair_gran_hooke.cpp
  rm -f ../pair_gran_hooke_history.cpp

  rm -f ../fix_freeze.h
  rm -f ../fix_pour.h
  rm -f ../fix_wall_gran.h
  rm -f ../pair_gran_hertz_history.h
  rm -f ../pair_gran_hooke.h
  rm -f ../pair_gran_hooke_history.h

fi
