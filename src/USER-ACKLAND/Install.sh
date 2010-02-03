# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_ackland_atom.cpp ..

  cp compute_ackland_atom.h ..

elif (test $1 = 0) then

  rm ../compute_ackland_atom.cpp

  rm ../compute_ackland_atom.h

fi
