# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

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

elif (test $1 = 0) then

  rm -f ../bond_class2.cpp
  rm -f ../angle_class2.cpp
  rm -f ../dihedral_class2.cpp
  rm -f ../improper_class2.cpp

  rm -f ../pair_lj_class2.cpp
  rm -f ../pair_lj_class2_coul_cut.cpp
  rm -f ../pair_lj_class2_coul_long.cpp

  rm -f ../bond_class2.h
  rm -f ../angle_class2.h
  rm -f ../dihedral_class2.h
  rm -f ../improper_class2.h

  rm -f ../pair_lj_class2.h
  rm -f ../pair_lj_class2_coul_cut.h
  rm -f ../pair_lj_class2_coul_long.h

fi
