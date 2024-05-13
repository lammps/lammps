# Install/Uninstall package files in LAMMPS
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

# package files without dependencies
action pair_sdpd_taitwater_isothermal.h
action pair_sdpd_taitwater_isothermal.cpp
action fix_meso_move.h
action fix_meso_move.cpp

# package files with dependencies
action fix_rigid_meso.h   fix_rigid.h
action fix_rigid_meso.cpp fix_rigid.h
