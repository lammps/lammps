# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp atom_vec_dipole.cpp ..
  cp pair_dipole_cut.cpp ..

  cp atom_vec_dipole.h ..
  cp pair_dipole_cut.h ..

elif (test $1 = 0) then

  rm ../atom_vec_dipole.cpp
  rm ../pair_dipole_cut.cpp

  rm ../atom_vec_dipole.h
  rm ../pair_dipole_cut.h

fi
