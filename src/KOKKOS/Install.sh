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

# force rebuild of files with LMP_KOKKOS switch

KOKKOS_INSTALLED=0
if (test -e ../Makefile.package) then
  KOKKOS_INSTALLED=`grep DLMP_KOKKOS ../Makefile.package | wc -l`
fi 

if (test $mode = 1) then
  if (test $KOKKOS_INSTALLED = 0) then
    touch ../accelerator_kokkos.h
  fi
elif (test $mode = 0) then
  if (test $KOKKOS_INSTALLED = 1) then
    touch ../accelerator_kokkos.h
  fi
fi

# list of files with optional dependcies

action angle_charmm_kokkos.cpp angle_charmm.cpp 
action angle_charmm_kokkos.h angle_charmm.h
action angle_class2_kokkos.cpp angle_class2.cpp 
action angle_class2_kokkos.h angle_class2.h
action angle_cosine_kokkos.cpp angle_cosine.cpp 
action angle_cosine_kokkos.h angle_cosine.h
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
action atom_vec_dpd_kokkos.cpp atom_vec_dpd.cpp
action atom_vec_dpd_kokkos.h atom_vec_dpd.h
action atom_vec_full_kokkos.cpp atom_vec_full.cpp
action atom_vec_full_kokkos.h atom_vec_full.h
action atom_vec_hybrid_kokkos.cpp
action atom_vec_hybrid_kokkos.h
action atom_vec_kokkos.cpp
action atom_vec_kokkos.h
action atom_vec_molecular_kokkos.cpp atom_vec_molecular.cpp
action atom_vec_molecular_kokkos.h atom_vec_molecular.h
action atom_vec_sphere_kokkos.cpp atom_vec_sphere.cpp
action atom_vec_sphere_kokkos.h atom_vec_sphere.h
action bond_class2_kokkos.cpp bond_class2.cpp 
action bond_class2_kokkos.h bond_class2.h
action bond_fene_kokkos.cpp bond_fene.cpp
action bond_fene_kokkos.h bond_fene.h
action bond_harmonic_kokkos.cpp bond_harmonic.cpp
action bond_harmonic_kokkos.h bond_harmonic.h
action comm_kokkos.cpp
action comm_kokkos.h
action comm_tiled_kokkos.cpp
action comm_tiled_kokkos.h
action compute_temp_kokkos.cpp
action compute_temp_kokkos.h
action dihedral_charmm_kokkos.cpp dihedral_charmm.cpp
action dihedral_charmm_kokkos.h dihedral_charmm.h
action dihedral_class2_kokkos.cpp dihedral_class2.cpp 
action dihedral_class2_kokkos.h dihedral_class2.h
action dihedral_harmonic_kokkos.cpp dihedral_harmonic.cpp 
action dihedral_harmonic_kokkos.h dihedral_harmonic.h
action dihedral_opls_kokkos.cpp dihedral_opls.cpp
action dihedral_opls_kokkos.h dihedral_opls.h
action domain_kokkos.cpp
action domain_kokkos.h
action fix_deform_kokkos.cpp
action fix_deform_kokkos.h
action fix_enforce2d_kokkos.cpp
action fix_enforce2d_kokkos.h
action fix_eos_table_rx_kokkos.cpp fix_eos_table_rx.cpp
action fix_eos_table_rx_kokkos.h fix_eos_table_rx.h  
action fix_freeze_kokkos.cpp fix_freeze.cpp
action fix_freeze_kokkos.h fix_freeze.h
action fix_gravity_kokkos.cpp
action fix_gravity_kokkos.h
action fix_langevin_kokkos.cpp
action fix_langevin_kokkos.h
action fix_minimize_kokkos.cpp
action fix_minimize_kokkos.h
action fix_neigh_history_kokkos.cpp
action fix_neigh_history_kokkos.h
action fix_nh_kokkos.cpp
action fix_nh_kokkos.h
action fix_nph_kokkos.cpp
action fix_nph_kokkos.h
action fix_npt_kokkos.cpp
action fix_npt_kokkos.h
action fix_nve_kokkos.cpp
action fix_nve_kokkos.h
action fix_nve_sphere_kokkos.cpp
action fix_nve_sphere_kokkos.h
action fix_nvt_kokkos.cpp
action fix_nvt_kokkos.h
action fix_property_atom_kokkos.cpp
action fix_property_atom_kokkos.h
action fix_qeq_reax_kokkos.cpp fix_qeq_reax.cpp
action fix_qeq_reax_kokkos.h fix_qeq_reax.h
action fix_reaxc_bonds_kokkos.cpp fix_reaxc_bonds.cpp
action fix_reaxc_bonds_kokkos.h fix_reaxc_bonds.h
action fix_reaxc_species_kokkos.cpp fix_reaxc_species.cpp
action fix_reaxc_species_kokkos.h fix_reaxc_species.h
action fix_rigid_kokkos.cpp fix_rigid.cpp
action fix_rigid_kokkos.h fix_rigid.h
action fix_setforce_kokkos.cpp
action fix_setforce_kokkos.h
action fix_shardlow_kokkos.cpp fix_shardlow.cpp
action fix_shardlow_kokkos.h fix_shardlow.h
action fix_momentum_kokkos.cpp
action fix_momentum_kokkos.h
action fix_wall_lj93_kokkos.cpp
action fix_wall_lj93_kokkos.h
action fix_wall_reflect_kokkos.cpp
action fix_wall_reflect_kokkos.h
action fix_dpd_energy_kokkos.cpp fix_dpd_energy.cpp
action fix_dpd_energy_kokkos.h fix_dpd_energy.h
action fix_rx_kokkos.cpp fix_rx.cpp
action fix_rx_kokkos.h fix_rx.h
action gridcomm_kokkos.cpp gridcomm.cpp
action gridcomm_kokkos.h gridcomm.h
action improper_class2_kokkos.cpp improper_class2.cpp 
action improper_class2_kokkos.h improper_class2.h
action improper_harmonic_kokkos.cpp improper_harmonic.cpp
action improper_harmonic_kokkos.h improper_harmonic.h
action kokkos.cpp
action kokkos.h
action kokkos_base.h
action kokkos_few.h
action kokkos_type.h
action memory_kokkos.h
action modify_kokkos.cpp
action modify_kokkos.h
action neigh_bond_kokkos.cpp
action neigh_bond_kokkos.h
action neigh_list_kokkos.cpp
action neigh_list_kokkos.h
action neighbor_kokkos.cpp
action neighbor_kokkos.h
action npair_copy_kokkos.cpp
action npair_copy_kokkos.h
action npair_halffull_kokkos.cpp
action npair_halffull_kokkos.h
action npair_skip_kokkos.cpp
action npair_skip_kokkos.h
action npair_kokkos.cpp
action npair_kokkos.h
action npair_ssa_kokkos.cpp npair_half_bin_newton_ssa.cpp
action npair_ssa_kokkos.h npair_half_bin_newton_ssa.h
action nbin_kokkos.cpp
action nbin_kokkos.h
action nbin_ssa_kokkos.cpp nbin_ssa.cpp
action nbin_ssa_kokkos.h nbin_ssa.h
action math_special_kokkos.cpp
action math_special_kokkos.h
action min_cg_kokkos.cpp
action min_cg_kokkos.h
action min_kokkos.cpp
action min_kokkos.h
action min_linesearch_kokkos.cpp
action min_linesearch_kokkos.h
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
action pair_dpd_fdt_energy_kokkos.cpp pair_dpd_fdt_energy.cpp
action pair_dpd_fdt_energy_kokkos.h pair_dpd_fdt_energy.h
action pair_eam_kokkos.cpp pair_eam.cpp
action pair_eam_kokkos.h pair_eam.h
action pair_eam_alloy_kokkos.cpp pair_eam_alloy.cpp
action pair_eam_alloy_kokkos.h pair_eam_alloy.h
action pair_eam_fs_kokkos.cpp pair_eam_fs.cpp
action pair_eam_fs_kokkos.h pair_eam_fs.h
action pair_exp6_rx_kokkos.cpp pair_exp6_rx.cpp
action pair_exp6_rx_kokkos.h pair_exp6_rx.h
action pair_gran_hooke_history_kokkos.h pair_gran_hooke_history.h
action pair_gran_hooke_history_kokkos.cpp pair_gran_hooke_history.cpp
action pair_hybrid_kokkos.cpp
action pair_hybrid_kokkos.h
action pair_hybrid_overlay_kokkos.cpp
action pair_hybrid_overlay_kokkos.h
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
action pair_morse_kokkos.cpp
action pair_morse_kokkos.h
action pair_multi_lucy_rx_kokkos.cpp pair_multi_lucy_rx.cpp
action pair_multi_lucy_rx_kokkos.h pair_multi_lucy_rx.h
action pair_reaxc_kokkos.cpp pair_reaxc.cpp
action pair_reaxc_kokkos.h pair_reaxc.h
action pair_snap_kokkos.cpp pair_snap.cpp
action pair_snap_kokkos.h pair_snap.h
action pair_snap_kokkos_impl.h pair_snap.cpp
action pair_sw_kokkos.cpp pair_sw.cpp
action pair_sw_kokkos.h pair_sw.h
action pair_vashishta_kokkos.cpp pair_vashishta.cpp
action pair_vashishta_kokkos.h pair_vashishta.h
action pair_table_kokkos.cpp
action pair_table_kokkos.h
action pair_table_rx_kokkos.cpp pair_table_rx.cpp
action pair_table_rx_kokkos.h pair_table_rx.h
action pair_tersoff_kokkos.cpp pair_tersoff.cpp
action pair_tersoff_kokkos.h pair_tersoff.h
action pair_tersoff_mod_kokkos.cpp pair_tersoff_mod.cpp
action pair_tersoff_mod_kokkos.h pair_tersoff_mod.h
action pair_tersoff_zbl_kokkos.cpp pair_tersoff_zbl.cpp
action pair_tersoff_zbl_kokkos.h pair_tersoff_zbl.h
action pair_yukawa_kokkos.cpp
action pair_yukawa_kokkos.h
action pair_zbl_kokkos.cpp
action pair_zbl_kokkos.h
action pppm_kokkos.cpp pppm.cpp
action pppm_kokkos.h pppm.h
action rand_pool_wrap_kokkos.cpp
action rand_pool_wrap_kokkos.h
action region_block_kokkos.cpp
action region_block_kokkos.h
action sna_kokkos.h sna.h
action sna_kokkos_impl.h sna.cpp
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
    sed -i -e '4 i \
CXX = $(CC)
' ../Makefile.package.settings
    sed -i -e '5 i \
include ..\/..\/lib\/kokkos\/Makefile.kokkos
' ../Makefile.package.settings
  fi

  #  comb/omp triggers a persistent bug in nvcc. deleting it.
  rm -f ../*_comb_omp.*

elif (test $1 = 2) then

  #  comb/omp triggers a persistent bug in nvcc. deleting it.
  rm -f ../*_comb_omp.*

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
