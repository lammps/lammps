# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp angle_cosine_shift.cpp ..
  cp angle_cosine_shift_exp.cpp ..
  cp bond_harmonic_shift.cpp ..
  cp bond_harmonic_shift_cut.cpp ..
  cp compute_temp_rotate.cpp ..
  cp dihedral_cosine_shift_exp.cpp ..
  cp fix_addtorque.cpp ..
  cp pair_dipole_sf.cpp ..

  cp angle_cosine_shift.h ..
  cp angle_cosine_shift_exp.h ..
  cp bond_harmonic_shift.h ..
  cp bond_harmonic_shift_cut.h ..
  cp compute_temp_rotate.h ..
  cp dihedral_cosine_shift_exp.h ..
  cp fix_addtorque.h ..
  cp pair_dipole_sf.h ..

elif (test $1 = 0) then

  rm -f ../angle_cosine_shift.cpp
  rm -f ../angle_cosine_shift_exp.cpp
  rm -f ../bond_harmonic_shift.cpp
  rm -f ../bond_harmonic_shift_cut.cpp
  rm -f ../compute_temp_rotate.cpp
  rm -f ../dihedral_cosine_shift_exp.cpp
  rm -f ../fix_addtorque.cpp
  rm -f ../pair_dipole_sf.cpp

  rm -f ../angle_cosine_shift.h
  rm -f ../angle_cosine_shift_exp.h
  rm -f ../bond_harmonic_shift.h
  rm -f ../bond_harmonic_shift_cut.h
  rm -f ../compute_temp_rotate.h
  rm -f ../dihedral_cosine_shift_exp.h
  rm -f ../fix_addtorque.h
  rm -f ../pair_dipole_sf.h

fi
