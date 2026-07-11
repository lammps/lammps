# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

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
action pair_dpd_coul_slater_long.cpp  pppm.cpp
action pair_dpd_coul_slater_long.h  pppm.h
action pair_dpd.cpp
action pair_dpd_ext.cpp
action pair_dpd_ext.h
action pair_dpd_ext_tstat.cpp
action pair_dpd_ext_tstat.h
action pair_dpd.h
action pair_dpd_tstat.cpp
action pair_dpd_tstat.h
