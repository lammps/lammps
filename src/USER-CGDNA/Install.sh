# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

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

# list of files with dependcies

action  bond_oxdna_fene.cpp bond_fene.h
action  bond_oxdna2_fene.cpp bond_fene.h
action  bond_oxdna_fene.h bond_fene.h
action  bond_oxdna2_fene.h bond_fene.h
action  fix_nve_dotc_langevin.cpp atom_vec_ellipsoid.h
action  fix_nve_dotc_langevin.h atom_vec_ellipsoid.h
action  fix_nve_dot.cpp atom_vec_ellipsoid.h
action  fix_nve_dot.h atom_vec_ellipsoid.h
action  mf_oxdna.h atom_vec_ellipsoid.h
action  pair_oxdna_coaxstk.cpp atom_vec_ellipsoid.h
action  pair_oxdna2_coaxstk.cpp atom_vec_ellipsoid.h
action  pair_oxdna_coaxstk.h atom_vec_ellipsoid.h
action  pair_oxdna2_coaxstk.h atom_vec_ellipsoid.h
action  pair_oxdna_excv.cpp atom_vec_ellipsoid.h
action  pair_oxdna2_excv.cpp atom_vec_ellipsoid.h
action  pair_oxdna_excv.h atom_vec_ellipsoid.h
action  pair_oxdna2_excv.h atom_vec_ellipsoid.h
action  pair_oxdna_hbond.cpp atom_vec_ellipsoid.h
action  pair_oxdna_hbond.h atom_vec_ellipsoid.h
action  pair_oxdna_stk.cpp atom_vec_ellipsoid.h
action  pair_oxdna2_stk.cpp atom_vec_ellipsoid.h
action  pair_oxdna_stk.h atom_vec_ellipsoid.h
action  pair_oxdna2_stk.h atom_vec_ellipsoid.h
action  pair_oxdna_xstk.cpp atom_vec_ellipsoid.h
action  pair_oxdna_xstk.h atom_vec_ellipsoid.h
action  pair_oxdna2_dh.cpp atom_vec_ellipsoid.h
action  pair_oxdna2_dh.h atom_vec_ellipsoid.h
