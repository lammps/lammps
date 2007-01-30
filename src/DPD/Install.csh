# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_dpd.h ..

  cp atom_vec_dpd.cpp ..
  cp pair_dpd.cpp ..

  cp atom_vec_dpd.h ..
  cp pair_dpd.h ..

else if ($1 == 0) then

  rm ../style_dpd.h
  touch ../style_dpd.h

  rm ../atom_vec_dpd.cpp
  rm ../pair_dpd.cpp

  rm ../atom_vec_dpd.h
  rm ../pair_dpd.h

endif
