# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_event_displace.cpp ..
  cp fix_event.cpp ..
  cp fix_neb.cpp ..
  cp neb.cpp ..
  cp prd.cpp ..
  cp temper.cpp ..

  cp compute_event_displace.h ..
  cp fix_event.h ..
  cp fix_neb.h ..
  cp neb.h ..
  cp prd.h ..
  cp temper.h ..

elif (test $1 = 0) then

  rm ../compute_event_displace.cpp
  rm ../fix_event.cpp
  rm ../fix_neb.cpp
  rm ../neb.cpp
  rm ../prd.cpp
  rm ../temper.cpp

  rm ../compute_event_displace.h
  rm ../fix_event.h
  rm ../fix_neb.h
  rm ../neb.h
  rm ../prd.h
  rm ../temper.h

fi
