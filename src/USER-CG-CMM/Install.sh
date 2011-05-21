# Install/unInstall package files in LAMMPS
# do not install child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../angle_harmonic.cpp) then
    cp angle_cg_cmm.h ..
    cp angle_cg_cmm.cpp ..
  fi

  if (test -e ../pppm.cpp) then
    cp pair_cg_cmm_coul_long.cpp ..
    cp pair_cg_cmm_coul_long.h ..
  fi

  cp cg_cmm_parms.h ..
  cp cg_cmm_parms.cpp ..

  cp pair_cmm_common.h ..
  cp pair_cmm_common.cpp ..
  cp pair_cg_cmm.cpp ..
  cp pair_cg_cmm.h ..
  cp pair_cg_cmm_coul_cut.cpp ..
  cp pair_cg_cmm_coul_cut.h ..

elif (test $1 = 0) then

  rm ../angle_cg_cmm.h
  rm ../angle_cg_cmm.cpp

  rm ../cg_cmm_parms.h
  rm ../cg_cmm_parms.cpp

  rm ../pair_cmm_common.h
  rm ../pair_cmm_common.cpp
  rm ../pair_cg_cmm.cpp
  rm ../pair_cg_cmm.h
  rm ../pair_cg_cmm_coul_cut.cpp
  rm ../pair_cg_cmm_coul_cut.h

  rm ../pair_cg_cmm_coul_long.cpp
  rm ../pair_cg_cmm_coul_long.h

fi
