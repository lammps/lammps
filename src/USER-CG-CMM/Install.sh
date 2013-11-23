# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# list of files with optional dependcies

action angle_sdk.cpp angle_harmonic.cpp
action angle_sdk.h angle_harmonic.cpp
action lj_sdk_common.h
action pair_lj_sdk.cpp
action pair_lj_sdk.h
action pair_lj_sdk_coul_long.cpp pppm.cpp
action pair_lj_sdk_coul_long.h pppm.cpp
action pair_lj_sdk_coul_msm.cpp msm.cpp
action pair_lj_sdk_coul_msm.h msm.cpp
