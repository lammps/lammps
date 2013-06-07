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

# list of files with optional dependencies

action atom_vec_dipole.h
action atom_vec_dipole.cpp
action pair_lj_cut_dipole_cut.h
action pair_lj_cut_dipole_cut.cpp

action pair_lj_cut_dipole_long.h ewald_disp.cpp
action pair_lj_cut_dipole_long.cpp ewald_disp.cpp
action pair_lj_long_dipole_long.h ewald_disp.cpp
action pair_lj_long_dipole_long.cpp ewald_disp.cpp
