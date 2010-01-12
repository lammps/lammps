# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp compute_ackland_atom.cpp ..

  cp compute_ackland_atom.h ..

else if ($1 == 0) then

  rm ../compute_ackland_atom.cpp

  rm ../compute_ackland_atom.h

endif
