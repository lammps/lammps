# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_wall_colloid.cpp ..
  cp pair_colloid.cpp ..
  cp pair_yukawa_colloid.cpp ..

  cp fix_wall_colloid.h ..
  cp pair_colloid.h ..
  cp pair_yukawa_colloid.h ..

elif (test $1 = 0) then

  rm -f ../fix_wall_colloid.cpp
  rm -f ../pair_colloid.cpp
  rm -f ../pair_yukawa_colloid.cpp

  rm -f ../fix_wall_colloid.h
  rm -f ../pair_colloid.h
  rm -f ../pair_yukawa_colloid.h

fi

