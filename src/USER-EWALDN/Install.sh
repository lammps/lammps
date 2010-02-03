# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p ewald_n.cpp ..
  cp -p pair_buck_coul.cpp ..
  cp -p pair_lj_coul.cpp ..

  cp -p ewald_n.h ..
  cp -p pair_buck_coul.h ..
  cp -p pair_lj_coul.h ..

  cp -p math_vector.h ..
  cp -p math_complex.h ..

elif (test $1 = 0) then

  rm ../ewald_n.cpp
  rm ../pair_buck_coul.cpp
  rm ../pair_lj_coul.cpp

  rm ../ewald_n.h
  rm ../pair_buck_coul.h
  rm ../pair_lj_coul.h

  rm ../math_vector.h
  rm ../math_complex.h

fi
