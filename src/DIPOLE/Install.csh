# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_dipole.h ..

  cp atom_vec_dipole.cpp ..
  cp compute_temp_dipole.cpp ..
  cp fix_nve_dipole.cpp ..
  cp pair_dipole_cut.cpp ..

  cp atom_vec_dipole.h ..
  cp compute_temp_dipole.h ..
  cp fix_nve_dipole.h ..
  cp pair_dipole_cut.h ..

else if ($1 == 0) then

  rm ../style_dipole.h
  touch ../style_dipole.h

  rm ../atom_vec_dipole.cpp
  rm ../compute_temp_dipole.cpp
  rm ../fix_nve_dipole.cpp
  rm ../pair_dipole_cut.cpp

  rm ../atom_vec_dipole.h
  rm ../compute_temp_dipole.h
  rm ../fix_nve_dipole.h
  rm ../pair_dipole_cut.h

endif
