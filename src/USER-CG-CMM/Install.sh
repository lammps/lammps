# Install/unInstall/Update package files in LAMMPS
# do not install or update child files if parent does not exist
mode=$1

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (test -z "$2" || test -e ../$2) then
    if (test $mode = 1) then
      cp $1 ..
    elif (test $mode = 2) then
      if (! cmp -s $1 ../$1) then
        echo "updating src/$1"
        cp $1 ..
      fi
    fi
  fi
}

#      package file               dependency (optional)
action angle_cg_cmm.h             angle_harmonic.cpp
action angle_cg_cmm.cpp           angle_harmonic.cpp

action angle_sdk.h                angle_harmonic.cpp
action angle_sdk.cpp              angle_harmonic.cpp

action pair_lj_sdk_coul_long.cpp  pppm.cpp
action pair_cg_cmm_coul_long.cpp  pppm.cpp
action pair_lj_sdk_coul_long.h    pppm.cpp
action pair_cg_cmm_coul_long.h    pppm.cpp

action pair_lj_sdk_coul_msm.cpp   msm.cpp
action pair_lj_sdk_coul_msm.h     msm.cpp

action cg_cmm_parms.h
action cg_cmm_parms.cpp

action pair_cmm_common.h
action pair_cmm_common.cpp
action pair_cg_cmm.cpp
action pair_cg_cmm.h
action pair_cg_cmm_coul_cut.cpp

action pair_lj_sdk.cpp
action pair_lj_sdk.h
action lj_sdk_common.h
