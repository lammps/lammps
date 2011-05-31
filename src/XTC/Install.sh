# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp dump_xtc.cpp ..

  cp dump_xtc.h ..

  cp xdr_compat.cpp ..
  cp xdr_compat.h ..

elif (test $1 = 0) then

  rm -f ../dump_xtc.cpp

  rm -f ../dump_xtc.h

  rm -f ../xdr_compat.cpp
  rm -f ../xdr_compat.h

fi
