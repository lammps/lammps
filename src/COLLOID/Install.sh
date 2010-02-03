# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp atom_vec_colloid.cpp ..
  cp fix_wall_colloid.cpp ..
  cp pair_colloid.cpp ..
  cp pair_lubricate.cpp ..
  cp pair_yukawa_colloid.cpp ..

  cp atom_vec_colloid.h ..
  cp fix_wall_colloid.h ..
  cp pair_colloid.h ..
  cp pair_lubricate.h ..
  cp pair_yukawa_colloid.h ..

elif (test $1 = 0) then

  rm ../atom_vec_colloid.cpp
  rm ../fix_wall_colloid.cpp
  rm ../pair_colloid.cpp
  rm ../pair_lubricate.cpp
  rm ../pair_yukawa_colloid.cpp

  rm ../atom_vec_colloid.h
  rm ../fix_wall_colloid.h
  rm ../pair_colloid.h
  rm ../pair_lubricate.h
  rm ../pair_yukawa_colloid.h

fi

