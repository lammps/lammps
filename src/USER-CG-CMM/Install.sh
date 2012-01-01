# Install/unInstall package files in LAMMPS
# do not install child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../angle_harmonic.cpp) then
    cp angle_cg_cmm.h ..
    cp angle_cg_cmm.cpp ..

    cp angle_sdk.h ..
    cp angle_sdk.cpp ..
  fi

  if (test -e ../pppm.cpp) then
    cp pair_cg_cmm_coul_long.cpp ..
    cp pair_cg_cmm_coul_long.h ..

    cp pair_lj_sdk_coul_long.cpp ..
    cp pair_lj_sdk_coul_long.h ..
  fi

  cp cg_cmm_parms.h ..
  cp cg_cmm_parms.cpp ..

  cp pair_cmm_common.h ..
  cp pair_cmm_common.cpp ..
  cp pair_cg_cmm.cpp ..
  cp pair_cg_cmm.h ..
  cp pair_cg_cmm_coul_cut.cpp ..
  cp pair_cg_cmm_coul_cut.h ..

  cp pair_lj_sdk.cpp ..
  cp pair_lj_sdk.h ..
  cp lj_sdk_common.h ..

elif (test $1 = 0) then

  rm -f ../angle_cg_cmm.h
  rm -f ../angle_cg_cmm.cpp

  rm -f ../cg_cmm_parms.h
  rm -f ../cg_cmm_parms.cpp

  rm -f ../pair_cmm_common.h
  rm -f ../pair_cmm_common.cpp
  rm -f ../pair_cg_cmm.cpp
  rm -f ../pair_cg_cmm.h
  rm -f ../pair_cg_cmm_coul_cut.cpp
  rm -f ../pair_cg_cmm_coul_cut.h

  rm -f ../pair_cg_cmm_coul_long.cpp
  rm -f ../pair_cg_cmm_coul_long.h

  rm -f ../lj_sdk_common.h

  rm -f ../angle_sdk.h
  rm -f ../angle_sdk.cpp

  rm -f ../pair_lj_sdk.cpp
  rm -f ../pair_lj_sdk.h
  rm -f ../pair_lj_sdk_coul_long.cpp
  rm -f ../pair_lj_sdk_coul_long.h
fi
