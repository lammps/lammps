# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_class2.h ..

  cp bond_class2.cpp ..
  cp angle_class2.cpp ..
  cp dihedral_class2.cpp ..
  cp improper_class2.cpp ..

  cp pair_lj_class2.cpp ..
  cp pair_lj_class2_coul_cut.cpp ..
  cp pair_lj_class2_coul_long.cpp ..

  cp bond_class2.h ..
  cp angle_class2.h ..
  cp dihedral_class2.h ..
  cp improper_class2.h ..

  cp pair_lj_class2.h ..
  cp pair_lj_class2_coul_cut.h ..
  cp pair_lj_class2_coul_long.h ..

else if ($1 == 0) then

  rm ../style_class2.h
  touch ../style_class2.h

  rm ../bond_class2.cpp
  rm ../angle_class2.cpp
  rm ../dihedral_class2.cpp
  rm ../improper_class2.cpp

  rm ../pair_lj_class2.cpp
  rm ../pair_lj_class2_coul_cut.cpp
  rm ../pair_lj_class2_coul_long.cpp

  rm ../bond_class2.h
  rm ../angle_class2.h
  rm ../dihedral_class2.h
  rm ../improper_class2.h

  rm ../pair_lj_class2.h
  rm ../pair_lj_class2_coul_cut.h
  rm ../pair_lj_class2_coul_long.h

endif
