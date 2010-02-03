# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_event_displace.cpp ..
  cp fix_event.cpp ..
  cp prd.cpp ..

  cp compute_event_displace.h ..
  cp fix_event.h ..
  cp prd.h ..

elif (test $1 = 0) then

  rm ../compute_event_displace.cpp
  rm ../fix_event.cpp
  rm ../prd.cpp

  rm ../compute_event_displace.h
  rm ../fix_event.h
  rm ../prd.h

fi
