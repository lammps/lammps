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

# list of files with optional dependencies

action angle_charmm_kokkos.cpp angle_charmm.cpp
action angle_charmm_kokkos.h angle_charmm.h
action angle_class2_kokkos.cpp angle_class2.cpp
action angle_class2_kokkos.h angle_class2.h
action angle_cosine_kokkos.cpp angle_cosine.cpp
action angle_cosine_kokkos.h angle_cosine.h
action angle_harmonic_kokkos.cpp angle_harmonic.cpp
action angle_harmonic_kokkos.h angle_harmonic.h
action angle_hybrid_kokkos.cpp angle_hybrid.cpp
action angle_hybrid_kokkos.h angle_hybrid.h
action angle_spica_kokkos.cpp angle_spica.cpp
action angle_spica_kokkos.h angle_spica.h
action atom_kokkos.cpp
action atom_kokkos.h
action atom_map_kokkos.cpp
action atom_vec_angle_kokkos.cpp atom_vec_angle.cpp
action atom_vec_angle_kokkos.h atom_vec_angle.h
action atom_vec_atomic_kokkos.cpp
action atom_vec_atomic_kokkos.h
action atom_vec_bond_kokkos.cpp atom_vec_bond.cpp
action atom_vec_bond_kokkos.h atom_vec_bond.h
action atom_vec_charge_kokkos.cpp
action atom_vec_charge_kokkos.h
action atom_vec_dipole_kokkos.cpp atom_vec_dipole.cpp
action atom_vec_dipole_kokkos.h atom_vec_dipole.h
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
action atom_vec_spin_kokkos.cpp atom_vec_spin.cpp
action atom_vec_spin_kokkos.h atom_vec_spin.h
action bond_class2_kokkos.cpp bond_class2.cpp
action bond_class2_kokkos.h bond_class2.h
action bond_fene_kokkos.cpp bond_fene.cpp
action bond_fene_kokkos.h bond_fene.h
action bond_harmonic_kokkos.cpp bond_harmonic.cpp
action bond_harmonic_kokkos.h bond_harmonic.h
action bond_hybrid_kokkos.cpp bond_hybrid.cpp
action bond_hybrid_kokkos.h bond_hybrid.h
action comm_brick_kokkos.cpp
action comm_brick_kokkos.h
action comm_brick_direct_kokkos.cpp
action comm_brick_direct_kokkos.h
action comm_tiled_kokkos.cpp
action comm_tiled_kokkos.h
action compute_ave_sphere_atom_kokkos.cpp compute_ave_sphere_atom.cpp
action compute_ave_sphere_atom_kokkos.h compute_ave_sphere_atom.h
action compute_coord_atom_kokkos.cpp
action compute_coord_atom_kokkos.h
action compute_erotate_sphere_kokkos.cpp
action compute_erotate_sphere_kokkos.h
action compute_composition_atom_kokkos.cpp compute_composition_atom.cpp
action compute_composition_atom_kokkos.h compute_composition_atom.h
action compute_orientorder_atom_kokkos.cpp
action compute_orientorder_atom_kokkos.h
action compute_temp_deform_kokkos.cpp
action compute_temp_deform_kokkos.h
action compute_temp_kokkos.cpp
action compute_temp_kokkos.h
action dihedral_charmm_kokkos.cpp dihedral_charmm.cpp
action dihedral_charmm_kokkos.h dihedral_charmm.h
action dihedral_charmmfsw_kokkos.cpp dihedral_charmmfsw.cpp
action dihedral_charmmfsw_kokkos.h dihedral_charmmfsw.h
action dihedral_class2_kokkos.cpp dihedral_class2.cpp
action dihedral_class2_kokkos.h dihedral_class2.h
action dihedral_harmonic_kokkos.cpp dihedral_harmonic.cpp
action dihedral_harmonic_kokkos.h dihedral_harmonic.h
action dihedral_opls_kokkos.cpp dihedral_opls.cpp
action dihedral_opls_kokkos.h dihedral_opls.h
action dihedral_hybrid_kokkos.cpp dihedral_hybrid.cpp
action dihedral_hybrid_kokkos.h dihedral_hybrid.h
action domain_kokkos.cpp
action domain_kokkos.h
action dynamical_matrix_kokkos.cpp dynamical_matrix.cpp
action dynamical_matrix_kokkos.h dynamical_matrix.h
action fft3d_kokkos.cpp fft3d.cpp
action fft3d_kokkos.h fft3d.h
action fftdata_kokkos.h fft3d.h
action fix_acks2_reaxff_kokkos.cpp fix_acks2_reaxff.cpp
action fix_acks2_reaxff_kokkos.h fix_acks2_reaxff.h
action fix_deform_kokkos.cpp
action fix_deform_kokkos.h
action fix_dpd_energy_kokkos.cpp fix_dpd_energy.cpp
action fix_dpd_energy_kokkos.h fix_dpd_energy.h
action fix_dt_reset_kokkos.cpp
action fix_dt_reset_kokkos.h
action fix_enforce2d_kokkos.cpp
action fix_enforce2d_kokkos.h
action fix_efield_kokkos.cpp
action fix_efield_kokkos.h
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
action fix_momentum_kokkos.cpp
action fix_momentum_kokkos.h
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
action fix_nvt_sllod_kokkos.cpp
action fix_nvt_sllod_kokkos.h
action fix_property_atom_kokkos.cpp
action fix_property_atom_kokkos.h
action fix_qeq_reaxff_kokkos.cpp fix_qeq_reaxff.cpp
action fix_qeq_reaxff_kokkos.h fix_qeq_reaxff.h
action fix_reaxff_bonds_kokkos.cpp fix_reaxff_bonds.cpp
action fix_reaxff_bonds_kokkos.h fix_reaxff_bonds.h
action compute_reaxff_atom_kokkos.cpp compute_reaxff_atom.cpp
action compute_reaxff_atom_kokkos.h compute_reaxff_atom.h
action fix_reaxff_species_kokkos.cpp fix_reaxff_species.cpp
action fix_reaxff_species_kokkos.h fix_reaxff_species.h
action fix_rx_kokkos.cpp fix_rx.cpp
action fix_rx_kokkos.h fix_rx.h
action fix_setforce_kokkos.cpp
action fix_setforce_kokkos.h
action fix_shake_kokkos.cpp fix_shake.cpp
action fix_shake_kokkos.h fix_shake.h
action fix_shardlow_kokkos.cpp fix_shardlow.cpp
action fix_shardlow_kokkos.h fix_shardlow.h
action fix_spring_self_kokkos.cpp
action fix_spring_self_kokkos.h
action fix_temp_berendsen_kokkos.cpp
action fix_temp_berendsen_kokkos.h
action fix_temp_rescale_kokkos.cpp
action fix_temp_rescale_kokkos.h
action fix_viscous_kokkos.cpp
action fix_viscous_kokkos.h
action fix_wall_flow_kokkos.cpp fix_wall_flow.cpp
action fix_wall_flow_kokkos.h fix_wall_flow.h
action fix_wall_gran_kokkos.cpp fix_wall_gran.cpp
action fix_wall_gran_kokkos.h fix_wall_gran.h
action fix_wall_gran_old.cpp fix_wall_gran.cpp
action fix_wall_gran_old.h fix_wall_gran.h
action fix_wall_lj93_kokkos.cpp
action fix_wall_lj93_kokkos.h
action fix_wall_reflect_kokkos.cpp
action fix_wall_reflect_kokkos.h
action grid3d_kokkos.cpp fft3d.h
action grid3d_kokkos.h fft3d.h
action improper_class2_kokkos.cpp improper_class2.cpp
action improper_class2_kokkos.h improper_class2.h
action improper_harmonic_kokkos.cpp improper_harmonic.cpp
action improper_harmonic_kokkos.h improper_harmonic.h
action improper_hybrid_kokkos.cpp improper_hybrid.cpp
action improper_hybrid_kokkos.h improper_hybrid.h
action kissfft_kokkos.h kissfft.h
action kokkos_base_fft.h fft3d.h
action kokkos_base.h
action kokkos_few.h
action kokkos_type.h
action kokkos.cpp
action kokkos.h
action math_special_kokkos.cpp
action math_special_kokkos.h
action meam_dens_final_kokkos.h meam_dens_final.cpp
action meam_dens_init_kokkos.h meam_dens_init.cpp
action meam_force_kokkos.h meam_force.cpp
action meam_funcs_kokkos.h meam_funcs.cpp
action meam_impl_kokkos.h meam_impl.cpp
action meam_kokkos.h meam.h
action meam_setup_done_kokkos.h meam_setup_done.cpp
action memory_kokkos.h
action min_cg_kokkos.cpp
action min_cg_kokkos.h
action min_kokkos.cpp
action min_kokkos.h
action min_linesearch_kokkos.cpp
action min_linesearch_kokkos.h
action mliap_data_kokkos.cpp mliap_data.cpp
action mliap_data_kokkos.h mliap_data.h
action mliap_descriptor_kokkos.h mliap_descriptor.h
action mliap_descriptor_so3_kokkos.cpp mliap_descriptor_so3.cpp
action mliap_descriptor_so3_kokkos.h mliap_descriptor_so3.h
action mliap_model_kokkos.h mliap_model.h
action mliap_model_linear_kokkos.cpp mliap_model_linear.cpp
action mliap_model_linear_kokkos.h mliap_model_linear.h
action mliap_model_python_kokkos.cpp mliap_model_linear.cpp
action mliap_model_python_kokkos.h mliap_model_linear.h
action mliap_so3_kokkos.cpp mliap_so3.cpp
action mliap_so3_kokkos.h mliap_so3.h
action mliap_unified_kokkos.cpp mliap_unified.cpp
action mliap_unified_kokkos.h mliap_unified.h
action modify_kokkos.cpp
action modify_kokkos.h
action nbin_kokkos.cpp
action nbin_kokkos.h
action nbin_ssa_kokkos.cpp nbin_ssa.cpp
action nbin_ssa_kokkos.h nbin_ssa.h
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
action npair_kokkos.cpp
action npair_kokkos.h
action npair_skip_kokkos.cpp
action npair_skip_kokkos.h
action npair_ssa_kokkos.cpp npair_half_bin_newton_ssa.cpp
action npair_ssa_kokkos.h npair_half_bin_newton_ssa.h
action npair_trim_kokkos.cpp
action npair_trim_kokkos.h
action pack_kokkos.h pack.h
action pair_adp_kokkos.cpp pair_adp.cpp
action pair_adp_kokkos.h pair_adp.h
action pair_brownian_kokkos.cpp pair_brownian.cpp
action pair_brownian_kokkos.h pair_brownian.h
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
action pair_dpd_ext_kokkos.cpp pair_dpd_ext.cpp
action pair_dpd_ext_kokkos.h pair_dpd_ext.h
action pair_dpd_ext_tstat_kokkos.cpp pair_dpd_ext_tstat.cpp
action pair_dpd_ext_tstat_kokkos.h pair_dpd_ext_tstat.h
action pair_dpd_fdt_energy_kokkos.cpp pair_dpd_fdt_energy.cpp
action pair_dpd_fdt_energy_kokkos.h pair_dpd_fdt_energy.h
action pair_dpd_kokkos.cpp pair_dpd.cpp
action pair_dpd_kokkos.h pair_dpd.h
action pair_dpd_tstat_kokkos.cpp pair_dpd_tstat.cpp
action pair_dpd_tstat_kokkos.h pair_dpd_tstat.h
action pair_eam_alloy_kokkos.cpp pair_eam_alloy.cpp
action pair_eam_alloy_kokkos.h pair_eam_alloy.h
action pair_eam_fs_kokkos.cpp pair_eam_fs.cpp
action pair_eam_fs_kokkos.h pair_eam_fs.h
action pair_eam_kokkos.cpp pair_eam.cpp
action pair_eam_kokkos.h pair_eam.h
action pair_exp6_rx_kokkos.cpp pair_exp6_rx.cpp
action pair_exp6_rx_kokkos.h pair_exp6_rx.h
action pair_gran_hooke_history_kokkos.cpp pair_gran_hooke_history.cpp
action pair_gran_hooke_history_kokkos.h pair_gran_hooke_history.h
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
action pair_lj_charmmfsw_coul_long_kokkos.cpp pair_lj_charmmfsw_coul_long.cpp
action pair_lj_charmmfsw_coul_long_kokkos.h pair_lj_charmmfsw_coul_long.h
action pair_lj_class2_coul_cut_kokkos.cpp pair_lj_class2_coul_cut.cpp
action pair_lj_class2_coul_cut_kokkos.h pair_lj_class2_coul_cut.h
action pair_lj_class2_coul_long_kokkos.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_kokkos.h pair_lj_class2_coul_long.h
action pair_lj_class2_kokkos.cpp pair_lj_class2.cpp
action pair_lj_class2_kokkos.h pair_lj_class2.h
action pair_lj_cut_coul_cut_kokkos.cpp
action pair_lj_cut_coul_cut_kokkos.h
action pair_lj_cut_coul_debye_kokkos.cpp pair_lj_cut_coul_debye.cpp
action pair_lj_cut_coul_debye_kokkos.h pair_lj_cut_coul_debye.h
action pair_lj_cut_coul_dsf_kokkos.cpp pair_lj_cut_coul_dsf.cpp
action pair_lj_cut_coul_dsf_kokkos.h pair_lj_cut_coul_dsf.h
action pair_lj_cut_coul_long_kokkos.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_kokkos.h pair_lj_cut_coul_long.h
action pair_lj_cut_dipole_cut_kokkos.cpp pair_lj_cut_dipole_cut.cpp
action pair_lj_cut_dipole_cut_kokkos.h pair_lj_cut_dipole_cut.h
action pair_lj_cut_kokkos.cpp
action pair_lj_cut_kokkos.h
action pair_lj_expand_coul_long_kokkos.cpp pair_lj_expand_coul_long.cpp
action pair_lj_expand_coul_long_kokkos.h pair_lj_expand_coul_long.h
action pair_lj_expand_kokkos.cpp
action pair_lj_expand_kokkos.h
action pair_lj_gromacs_coul_gromacs_kokkos.cpp pair_lj_gromacs_coul_gromacs.cpp
action pair_lj_gromacs_coul_gromacs_kokkos.h pair_lj_gromacs_coul_gromacs.h
action pair_lj_gromacs_kokkos.cpp pair_lj_gromacs.cpp
action pair_lj_gromacs_kokkos.h pair_lj_gromacs.h
action pair_lj_spica_coul_long_kokkos.cpp pair_lj_spica_coul_long.cpp
action pair_lj_spica_coul_long_kokkos.h pair_lj_spica_coul_long.h
action pair_lj_spica_kokkos.cpp pair_lj_spica.cpp
action pair_lj_spica_kokkos.h pair_lj_spica.h
action pair_meam_kokkos.cpp pair_meam.cpp
action pair_meam_kokkos.h pair_meam.h
action pair_meam_ms_kokkos.cpp pair_meam_ms.cpp
action pair_meam_ms_kokkos.h pair_meam_ms.h
action pair_mliap_kokkos.cpp pair_mliap.cpp
action pair_mliap_kokkos.h pair_mliap.h
action pair_morse_kokkos.cpp
action pair_morse_kokkos.h
action pair_multi_lucy_rx_kokkos.cpp pair_multi_lucy_rx.cpp
action pair_multi_lucy_rx_kokkos.h pair_multi_lucy_rx.h
action pair_pace_extrapolation_kokkos.cpp pair_pace_extrapolation.cpp
action pair_pace_extrapolation_kokkos.h pair_pace_extrapolation.h
action pair_pod_kokkos.cpp pair_pod.cpp
action pair_pod_kokkos.h pair_pod.h
action pair_pace_kokkos.cpp pair_pace.cpp
action pair_pace_kokkos.h pair_pace.h
action pair_reaxff_kokkos.cpp pair_reaxff.cpp
action pair_reaxff_kokkos.h pair_reaxff.h
action pair_snap_kokkos_impl.h pair_snap.cpp
action pair_snap_kokkos.cpp pair_snap.cpp
action pair_snap_kokkos.h pair_snap.h
action pair_soft_kokkos.cpp
action pair_soft_kokkos.h
action pair_sw_kokkos.cpp pair_sw.cpp
action pair_sw_kokkos.h pair_sw.h
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
action pair_uf3_kokkos.cpp pair_uf3.cpp
action pair_uf3_kokkos.h pair_uf3.h
action pair_vashishta_kokkos.cpp pair_vashishta.cpp
action pair_vashishta_kokkos.h pair_vashishta.h
action pair_yukawa_kokkos.cpp
action pair_yukawa_kokkos.h
action pair_yukawa_colloid_kokkos.cpp pair_yukawa_colloid.cpp
action pair_yukawa_colloid_kokkos.h pair_yukawa_colloid.h
action pair_zbl_kokkos.cpp
action pair_zbl_kokkos.h
action pppm_kokkos.cpp pppm.cpp
action pppm_kokkos.h pppm.h
action rand_pool_wrap_kokkos.cpp
action rand_pool_wrap_kokkos.h
action region_block_kokkos.cpp
action region_block_kokkos.h
action remap_kokkos.cpp remap.cpp
action remap_kokkos.h remap.h
action sna_kokkos_impl.h sna.cpp
action sna_kokkos.h sna.h
action third_order_kokkos.cpp dynamical_matrix.cpp
action third_order_kokkos.h dynamical_matrix.h
action transpose_helper_kokkos.h
action verlet_kokkos.cpp
action verlet_kokkos.h

# Install cython pyx file only if non-KOKKOS version is present
action mliap_model_python_couple_kokkos.pyx mliap_model_python_couple.pyx
action mliap_unified_couple_kokkos.pyx mliap_unified_couple.pyx

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
    sed -i -e '/^[ \t]*include.*kokkos.*$/d' ../Makefile.package.settings
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
    sed -i -e '/^[ \t]*include.*kokkos.*$/d' ../Makefile.package.settings
  fi

fi

# Python cython stuff. Only need to convert/remove sources.
# Package settings were already done in ML-IAP package Install.sh script.

if (test $1 = 1) then
  if (type cythonize > /dev/null 2>&1 && test -e ../python_impl.cpp) then
    cythonize -3 ../mliap_model_python_couple_kokkos.pyx
    cythonize -3 ../mliap_unified_couple_kokkos.pyx
  fi

elif (test $1 = 0) then
  rm -f ../mliap_model_python_couple_kokkos.cpp ../mliap_model_python_couple_kokkos.h
  rm -f ../mliap_unified_couple_kokkos.cpp ../mliap_unified_couple_kokkos.h

elif (test $1 = 2) then
  if (type cythonize > /dev/null 2>&1 && test -e ../python_impl.cpp) then
    cythonize -3 ../mliap_model_python_couple_kokkos.pyx
    cythonize -3 ../mliap_unified_couple_kokkos.pyx
  else
    rm -f ../mliap_model_python_couple_kokkos.cpp ../mliap_model_python_couple_kokkos.h
    rm -f ../mliap_unified_couple_kokkos.cpp ../mliap_unified_couple_kokkos.h
  fi
fi
