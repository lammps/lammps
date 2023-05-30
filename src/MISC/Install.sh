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
action bond_special.cpp
action bond_special.h
action compute_viscosity_cos.cpp
action compute_viscosity_cos.h
action fix_accelerate_cos.cpp
action fix_accelerate_cos.h
action fix_imd.cpp
action fix_imd.h
action fix_ipi.cpp
action fix_ipi.h
action fix_srp.cpp
action fix_srp.h
action pair_agni.cpp
action pair_agni.h
action pair_list.cpp
action pair_list.h
action pair_srp.cpp
action pair_srp.h
action pair_tracker.cpp
action pair_tracker.h

# package files with dependencies
action pair_srp_react.cpp fix_bond_break.h
action pair_srp_react.h fix_bond_break.h
action fix_srp_react.cpp fix_bond_break.h
action fix_srp_react.h fix_bond_break.h
