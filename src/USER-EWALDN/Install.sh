# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p ewald_n.cpp ..
  cp -p pair_buck_coul.cpp ..
  #cp -p pair_lj_dipole.cpp ..
  cp -p pair_lj_coul.cpp ..

  cp -p ewald_n.h ..
  cp -p pair_buck_coul.h ..
  #cp -p pair_lj_dipole.h ..
  cp -p pair_lj_coul.h ..

  cp -p math_vector.h ..
  cp -p math_complex.h ..

elif (test $1 = 0) then

  rm -f ../ewald_n.cpp
  rm -f ../pair_buck_coul.cpp
  #rm -f ../pair_lj_dipole.cpp
  rm -f ../pair_lj_coul.cpp

  rm -f ../ewald_n.h
  rm -f ../pair_buck_coul.h
  #rm -f ../pair_lj_dipole.h
  rm -f ../pair_lj_coul.h

  rm -f ../math_vector.h
  rm -f ../math_complex.h

fi
