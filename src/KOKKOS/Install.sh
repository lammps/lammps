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

action angle_charmm_kokkos.cpp angle_charmm.cpp 
action angle_charmm_kokkos.h angle_charmm.h
action angle_harmonic_kokkos.cpp angle_harmonic.cpp 
action angle_harmonic_kokkos.h angle_harmonic.h 
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
action bond_fene_kokkos.cpp bond_fene.cpp
action bond_fene_kokkos.h bond_fene.h
action bond_harmonic_kokkos.cpp bond_harmonic.cpp
action bond_harmonic_kokkos.h bond_harmonic.h
action comm_kokkos.cpp
action comm_kokkos.h
action compute_temp_kokkos.cpp
action compute_temp_kokkos.h
action dihedral_charmm_kokkos.cpp dihedral_charmm.cpp
action dihedral_charmm_kokkos.h dihedral_charmm.h
action dihedral_opls_kokkos.cpp dihedral_opls.cpp
action dihedral_opls_kokkos.h dihedral_opls.h
action domain_kokkos.cpp
action domain_kokkos.h
action fix_deform_kokkos.cpp
action fix_deform_kokkos.h
action fix_langevin_kokkos.cpp
action fix_langevin_kokkos.h
action fix_nh_kokkos.cpp
action fix_nh_kokkos.h
action fix_nph_kokkos.cpp
action fix_nph_kokkos.h
action fix_npt_kokkos.cpp
action fix_npt_kokkos.h
action fix_nve_kokkos.cpp
action fix_nve_kokkos.h
action fix_nvt_kokkos.cpp
action fix_nvt_kokkos.h
action fix_setforce_kokkos.cpp
action fix_setforce_kokkos.h
action fix_wall_reflect_kokkos.cpp
action fix_wall_reflect_kokkos.h
action improper_harmonic_kokkos.cpp improper_harmonic.cpp
action improper_harmonic_kokkos.h improper_harmonic.h
action kokkos.cpp
action kokkos.h
action kokkos_type.h
action memory_kokkos.h
action modify_kokkos.cpp
action modify_kokkos.h
action neigh_bond_kokkos.cpp
action neigh_bond_kokkos.h
action neigh_full_kokkos.h
action neigh_list_kokkos.cpp
action neigh_list_kokkos.h
action neighbor_kokkos.cpp
action neighbor_kokkos.h
action pair_buck_coul_cut_kokkos.cpp
action pair_buck_coul_cut_kokkos.h
action pair_buck_coul_long_kokkos.cpp pair_buck_coul_long.cpp
action pair_buck_coul_long_kokkos.h pair_buck_coul_long.h
action pair_buck_kokkos.cpp
action pair_buck_kokkos.h
action pair_coul_cut_kokkos.cpp
action pair_coul_cut_kokkos.h
action pair_coul_debye_kokkos.cpp
action pair_coul_debye_kokkos.h
action pair_coul_dsf_kokkos.cpp
action pair_coul_dsf_kokkos.h
action pair_coul_long_kokkos.cpp pair_coul_long.cpp
action pair_coul_long_kokkos.h pair_coul_long.h
action pair_coul_wolf_kokkos.cpp
action pair_coul_wolf_kokkos.h
action pair_eam_kokkos.cpp pair_eam.cpp
action pair_eam_kokkos.h pair_eam.h
action pair_eam_alloy_kokkos.cpp pair_eam_alloy.cpp
action pair_eam_alloy_kokkos.h pair_eam_alloy.h
action pair_eam_fs_kokkos.cpp pair_eam_fs.cpp
action pair_eam_fs_kokkos.h pair_eam_fs.h
action pair_kokkos.h
action pair_lj_charmm_coul_charmm_implicit_kokkos.cpp pair_lj_charmm_coul_charmm_implicit.cpp
action pair_lj_charmm_coul_charmm_implicit_kokkos.h pair_lj_charmm_coul_charmm_implicit.h
action pair_lj_charmm_coul_charmm_kokkos.cpp pair_lj_charmm_coul_charmm.cpp
action pair_lj_charmm_coul_charmm_kokkos.h pair_lj_charmm_coul_charmm.h
action pair_lj_charmm_coul_long_kokkos.cpp pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_long_kokkos.h pair_lj_charmm_coul_long.h
action pair_lj_class2_coul_cut_kokkos.cpp pair_lj_class2_coul_cut.cpp
action pair_lj_class2_coul_cut_kokkos.h pair_lj_class2_coul_cut.h
action pair_lj_class2_coul_long_kokkos.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_kokkos.h pair_lj_class2_coul_long.h
action pair_lj_class2_kokkos.cpp pair_lj_class2.cpp
action pair_lj_class2_kokkos.h pair_lj_class2.h
action pair_lj_cut_coul_cut_kokkos.cpp
action pair_lj_cut_coul_cut_kokkos.h
action pair_lj_cut_coul_debye_kokkos.cpp
action pair_lj_cut_coul_debye_kokkos.h
action pair_lj_cut_coul_dsf_kokkos.cpp
action pair_lj_cut_coul_dsf_kokkos.h
action pair_lj_cut_coul_long_kokkos.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_kokkos.h pair_lj_cut_coul_long.h
action pair_lj_cut_kokkos.cpp
action pair_lj_cut_kokkos.h
action pair_lj_expand_kokkos.cpp
action pair_lj_expand_kokkos.h
action pair_lj_gromacs_coul_gromacs_kokkos.cpp
action pair_lj_gromacs_coul_gromacs_kokkos.h
action pair_lj_gromacs_kokkos.cpp
action pair_lj_gromacs_kokkos.h
action pair_lj_sdk_kokkos.cpp pair_lj_sdk.cpp
action pair_lj_sdk_kokkos.h pair_lj_sdk.h
action pair_sw_kokkos.cpp pair_sw.cpp
action pair_sw_kokkos.h pair_sw.h
action pair_table_kokkos.cpp
action pair_table_kokkos.h
action pair_tersoff_kokkos.cpp pair_tersoff.cpp
action pair_tersoff_kokkos.h pair_tersoff.h
action pair_tersoff_mod_kokkos.cpp pair_tersoff_mod.cpp
action pair_tersoff_mod_kokkos.h pair_tersoff_mod.h
action pair_tersoff_zbl_kokkos.cpp pair_tersoff_zbl.cpp
action pair_tersoff_zbl_kokkos.h pair_tersoff_zbl.h
action region_block_kokkos.cpp
action region_block_kokkos.h
action verlet_kokkos.cpp
action verlet_kokkos.h

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_KOKKOS |' ../Makefile.package
#    sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/kokkos\/core\/src |' ../Makefile.package
    sed -i -e 's|^PKG_CPP_DEPENDS =[ \t]*|&$(KOKKOS_CPP_DEPENDS) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&$(KOKKOS_LIBS) |' ../Makefile.package
    sed -i -e 's|^PKG_LINK_DEPENDS =[ \t]*|&$(KOKKOS_LINK_DEPENDS) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(KOKKOS_LDFLAGS) |' ../Makefile.package
#    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(kokkos_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/CXX\ =\ \$(CC)/d' ../Makefile.package.settings
    sed -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \CXX = $(CC)' ../Makefile.package.settings
    sed -i -e '5 i \include ..\/..\/lib\/kokkos\/Makefile.kokkos' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kokkos[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*KOKKOS[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/CXX\ =\ \$(CC)/d' ../Makefile.package.settings
    sed -i -e '/^include.*kokkos.*$/d' ../Makefile.package.settings
  fi

fi
