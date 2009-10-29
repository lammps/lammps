# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_prd.h ..

  cp compute_event_displace.cpp ..
  cp fix_event.cpp ..
  cp prd.cpp ..

  cp compute_event_displace.h ..
  cp fix_event.h ..
  cp prd.h ..

else if ($1 == 0) then

  rm ../style_prd.h
  touch ../style_prd.h

  rm ../compute_event_displace.cpp
  rm ../fix_event.cpp
  rm ../prd.cpp

  rm ../compute_event_displace.h
  rm ../fix_event.h
  rm ../prd.h

endif
