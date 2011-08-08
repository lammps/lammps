# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_temp_rotate.cpp ..
  cp fix_addtorque.cpp ..

  cp compute_temp_rotate.h ..
  cp fix_addtorque.h ..

elif (test $1 = 0) then

  rm -f ../compute_temp_rotate.cpp
  rm -f ../fix_addtorque.cpp

  rm -f ../compute_temp_rotate.h
  rm -f ../fix_addtorque.h

fi
