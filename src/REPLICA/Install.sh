# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_event_displace.cpp ..
  cp fix_event.cpp ..
  cp fix_event_prd.cpp ..
  cp fix_event_tad.cpp ..
  cp fix_neb.cpp ..
  cp neb.cpp ..
  cp prd.cpp ..
  cp tad.cpp ..
  cp temper.cpp ..

  cp compute_event_displace.h ..
  cp fix_event.h ..
  cp fix_event_prd.h ..
  cp fix_event_tad.h ..
  cp fix_neb.h ..
  cp neb.h ..
  cp prd.h ..
  cp tad.h ..
  cp temper.h ..

elif (test $1 = 0) then

  rm ../compute_event_displace.cpp
  rm ../fix_event.cpp
  rm ../fix_event_prd.cpp
  rm ../fix_event_tad.cpp
  rm ../fix_neb.cpp
  rm ../neb.cpp
  rm ../prd.cpp
  rm ../tad.cpp
  rm ../temper.cpp

  rm ../compute_event_displace.h
  rm ../fix_event.h
  rm ../fix_event_prd.h
  rm ../fix_event_tad.h
  rm ../fix_neb.h
  rm ../neb.h
  rm ../prd.h
  rm ../tad.h
  rm ../temper.h

fi
