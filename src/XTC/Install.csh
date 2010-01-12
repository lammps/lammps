# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp dump_xtc.cpp ..

  cp dump_xtc.h ..

  cp xdr_compat.cpp ..
  cp xdr_compat.h ..

else if ($1 == 0) then

  rm ../dump_xtc.cpp

  rm ../dump_xtc.h

  rm ../xdr_compat.cpp
  rm ../xdr_compat.h

endif
