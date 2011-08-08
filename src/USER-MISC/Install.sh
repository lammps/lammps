# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp angle_cosineshift.cpp ..
  cp angle_cosineshiftexp.cpp ..
  cp bond_harmonic_shift.cpp ..
  cp bond_harmonic_shift_cut.cpp ..
  cp compute_temp_rotate.cpp ..
  cp dihedral_cosineshiftexp.cpp ..
  cp fix_addtorque.cpp ..
  cp pair_dipole_sf.cpp ..

  cp angle_cosineshift.h ..
  cp angle_cosineshiftexp.h ..
  cp bond_harmonic_shift.h ..
  cp bond_harmonic_shift_cut.h ..
  cp compute_temp_rotate.h ..
  cp dihedral_cosineshiftexp.h ..
  cp fix_addtorque.h ..
  cp pair_dipole_sf.h ..

elif (test $1 = 0) then

  rm -f ../angle_cosineshift.cpp
  rm -f ../angle_cosineshiftexp.cpp
  rm -f ../bond_harmonic_shift.cpp
  rm -f ../bond_harmonic_shift_cut.cpp
  rm -f ../compute_temp_rotate.cpp
  rm -f ../dihedral_cosineshiftexp.cpp
  rm -f ../fix_addtorque.cpp
  rm -f ../pair_dipole_sf.cpp

  rm -f ../angle_cosineshift.h
  rm -f ../angle_cosineshiftexp.h
  rm -f ../bond_harmonic_shift.h
  rm -f ../bond_harmonic_shift_cut.h
  rm -f ../compute_temp_rotate.h
  rm -f ../dihedral_cosineshiftexp.h
  rm -f ../fix_addtorque.h
  rm -f ../pair_dipole_sf.h

fi
