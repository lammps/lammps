# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp compute_event_displace.cpp ..
  cp fix_event.cpp ..
  cp fix_event_prd.cpp ..
  cp fix_event_tad.cpp ..
  cp fix_neb.cpp ..
  cp verlet_split.cpp ..
  cp neb.cpp ..
  cp prd.cpp ..
  cp tad.cpp ..
  cp temper.cpp ..

  cp compute_event_displace.h ..
  cp fix_event.h ..
  cp fix_event_prd.h ..
  cp fix_event_tad.h ..
  cp fix_neb.h ..
  cp verlet_split.h ..
  cp neb.h ..
  cp prd.h ..
  cp tad.h ..
  cp temper.h ..

elif (test $1 = 0) then

  rm -f ../compute_event_displace.cpp
  rm -f ../fix_event.cpp
  rm -f ../fix_event_prd.cpp
  rm -f ../fix_event_tad.cpp
  rm -f ../fix_neb.cpp
  rm -f ../verlet_split.cpp
  rm -f ../neb.cpp
  rm -f ../prd.cpp
  rm -f ../tad.cpp
  rm -f ../temper.cpp

  rm -f ../compute_event_displace.h
  rm -f ../fix_event.h
  rm -f ../fix_event_prd.h
  rm -f ../fix_event_tad.h
  rm -f ../fix_neb.h
  rm -f ../verlet_split.h
  rm -f ../neb.h
  rm -f ../prd.h
  rm -f ../tad.h
  rm -f ../temper.h

fi
