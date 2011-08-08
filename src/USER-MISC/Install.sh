# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_temp_rotate.cpp ..
  cp fix_addtorque.cpp ..
  cp pair_dipole_sf.cpp ..

  cp compute_temp_rotate.h ..
  cp fix_addtorque.h ..
  cp pair_dipole_sf.h ..

elif (test $1 = 0) then

  rm -f ../compute_temp_rotate.cpp
  rm -f ../fix_addtorque.cpp
  rm -f ../pair_dipole_sf.cpp

  rm -f ../compute_temp_rotate.h
  rm -f ../fix_addtorque.h
  rm -f ../pair_dipole_sf.h

fi
