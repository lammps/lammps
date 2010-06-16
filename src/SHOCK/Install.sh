# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_msst.cpp ..
  cp compute_vsum.cpp ..

  cp fix_msst.h ..
  cp compute_vsum.h ..

elif (test $1 = 0) then

  rm ../fix_msst.cpp
  rm ../compute_vsum.cpp

  rm ../fix_msst.h
  rm ../compute_vsum.h

fi
