# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp -p style_user_ewaldn.h ..

  cp -p ewald_n.cpp ..
  cp -p pair_buck_coul.cpp ..
  cp -p pair_lj_coul.cpp ..

  cp -p ewald_n.h ..
  cp -p pair_buck_coul.h ..
  cp -p pair_lj_coul.h ..

  cp -p math_vector.h ..
  cp -p math_complex.h ..

else if ($1 == 0) then

  rm ../style_user_ewaldn.h
  touch ../style_user_ewaldn.h

  rm ../ewald_n.cpp
  rm ../pair_buck_coul.cpp
  rm ../pair_lj_coul.cpp

  rm ../ewald_n.h
  rm ../pair_buck_coul.h
  rm ../pair_lj_coul.h

  rm ../math_vector.h
  rm ../math_complex.h

endif
