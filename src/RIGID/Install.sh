# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp fix_rigid.cpp ..
  cp fix_rigid_nh.cpp ..
  cp fix_rigid_nph.cpp ..
  cp fix_rigid_npt.cpp ..
  cp fix_rigid_nve.cpp ..
  cp fix_rigid_nvt.cpp ..
  cp fix_rigid_small.cpp ..

  cp fix_rigid.h ..
  cp fix_rigid_nh.h ..
  cp fix_rigid_nph.h ..
  cp fix_rigid_npt.h ..
  cp fix_rigid_nve.h ..
  cp fix_rigid_nvt.h ..
  cp fix_rigid_small.h ..

elif (test $1 = 0) then

  rm -f ../fix_rigid.cpp
  rm -f ../fix_rigid_nh.cpp
  rm -f ../fix_rigid_nph.cpp
  rm -f ../fix_rigid_npt.cpp
  rm -f ../fix_rigid_nve.cpp
  rm -f ../fix_rigid_nvt.cpp
  rm -f ../fix_rigid_small.cpp

  rm -f ../fix_rigid.h
  rm -f ../fix_rigid_nh.h
  rm -f ../fix_rigid_nph.h
  rm -f ../fix_rigid_npt.h
  rm -f ../fix_rigid_nve.h
  rm -f ../fix_rigid_nvt.h
  rm -f ../fix_rigid_small.h

fi

