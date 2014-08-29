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

# force rebuild of files with LMP_KOKKOS switch

touch ../accelerator_kokkos.h
touch ../memory.h

# list of files with optional dependcies

action atom_kokkos.cpp
action atom_kokkos.h
action atom_vec_angle_kokkos.cpp atom_vec_angle.cpp
action atom_vec_angle_kokkos.h atom_vec_angle.h
action atom_vec_atomic_kokkos.cpp
action atom_vec_atomic_kokkos.h
action atom_vec_bond_kokkos.cpp atom_vec_bond.cpp
action atom_vec_bond_kokkos.h atom_vec_bond.h
action atom_vec_charge_kokkos.cpp
action atom_vec_charge_kokkos.h
action atom_vec_full_kokkos.cpp atom_vec_full.cpp
action atom_vec_full_kokkos.h atom_vec_full.h
action atom_vec_kokkos.cpp
action atom_vec_kokkos.h
action atom_vec_molecular_kokkos.cpp atom_vec_molecular.cpp
action atom_vec_molecular_kokkos.h atom_vec_molecular.h
action comm_kokkos.cpp
action comm_kokkos.h
action domain_kokkos.cpp
action domain_kokkos.h
action fix_langevin_kokkos.cpp
action fix_langevin_kokkos.h
action fix_nve_kokkos.cpp
action fix_nve_kokkos.h
action kokkos.cpp
action kokkos.h
action kokkos_type.h
action memory_kokkos.h
action modify_kokkos.cpp
action modify_kokkos.h
action neigh_full_kokkos.h
action neigh_list_kokkos.cpp
action neigh_list_kokkos.h
action neighbor_kokkos.cpp
action neighbor_kokkos.h
action pair_coul_cut_kokkos.cpp
action pair_coul_cut_kokkos.h
action pair_kokkos.h
action pair_lj_cut_coul_cut_kokkos.cpp
action pair_lj_cut_coul_cut_kokkos.h
action pair_lj_cut_coul_long_kokkos.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_kokkos.h pair_lj_cut_coul_long.h
action pair_lj_cut_kokkos.cpp
action pair_lj_cut_kokkos.h
action pair_table_kokkos.cpp
action pair_table_kokkos.h
action verlet_kokkos.cpp
action verlet_kokkos.h

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_KOKKOS |' ../Makefile.package
#    sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/kokkos\/core\/src |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lkokkoscore |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(KOKKOS_INC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(KOKKOS_LINK) |' ../Makefile.package
#    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(kokkos_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/kokkos\/Makefile.lammps
' ../Makefile.package.settings

  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
  fi

fi
