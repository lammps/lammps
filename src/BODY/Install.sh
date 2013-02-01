# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp body_nparticle.cpp ..
  cp compute_body_local.cpp ..
  cp fix_nve_body.cpp ..
  cp pair_body.cpp ..

  cp body_nparticle.h ..
  cp compute_body_local.h ..
  cp fix_nve_body.h ..
  cp pair_body.h ..

elif (test $1 = 0) then

  rm -f ../body_nparticle.cpp
  rm -f ../compute_body_local.cpp
  rm -f ../fix_nve_body.cpp
  rm -f ../pair_body.cpp

  rm -f ../body_nparticle.h
  rm -f ../compute_body_local.h
  rm -f ../fix_nve_body.h
  rm -f ../pair_body.h

fi
