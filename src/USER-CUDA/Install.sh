# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude CUDA library
# do not copy child files if parent does not exist

if (test $1 = 1) then

  if (test ! -e ../Makefile.package) then
    cp ../Makefile.package.empty ../Makefile.package
  fi

  sed -i -e '/^include.*cuda.*$/d' ../Makefile.package
  sed -i -e 's/[^ \t]*cuda[^ \t]* //g' ../Makefile.package
  sed -i -e 's/[^ \t]*CUDA[^ \t]* //g' ../Makefile.package
  sed -i -e 's/[^ \t]*lrt[^ \t]* //g' ../Makefile.package
  sed -i '4 i include ..\/..\/lib\/cuda\/Makefile.common' ../Makefile.package
  sed -i -e 's|^PKG_INC =[ \t]*|&-I..\/..\/lib\/cuda -DLMP_USER_CUDA |' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/cuda |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-llammpscuda |' ../Makefile.package
  sed -i -e 's|^PKG_SYSINC =[ \t]*|&-I$(CUDA_INSTALL_PATH)\/include |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&-L$(CUDA_INSTALL_PATH)\/lib64 -L$(CUDA_INSTALL_PATH)\/lib $(CUDA_USRLIB_CONDITIONAL) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&-lcuda -lcudart -lrt |' ../Makefile.package
 
  # force rebuild of dependencies since adding -DLMP_USER_CUDA switch

  touch ../accelerator_cuda.h

  if (test -e ../atom_vec_angle.cpp) then
    cp atom_vec_angle_cuda.cpp ..
    cp atom_vec_angle_cuda.h ..
  fi

  if (test -e ../atom_vec_full.cpp) then
    cp atom_vec_full_cuda.cpp ..
    cp atom_vec_full_cuda.h ..
  fi

  if (test -e ../fix_freeze.cpp) then
    cp fix_freeze_cuda.cpp ..
    cp fix_freeze_cuda.h ..
  fi

  if (test -e ../pair_born_coul_long.cpp) then
    cp pair_born_coul_long_cuda.cpp ..
    cp pair_born_coul_long_cuda.h ..
  fi

  if (test -e ../pair_buck_coul_long.cpp) then
    cp pair_buck_coul_long_cuda.cpp ..
    cp pair_buck_coul_long_cuda.h ..
  fi

  if (test -e ../pair_cg_cmm.cpp) then
    cp pair_cg_cmm_cuda.cpp ..
    cp pair_cg_cmm_coul_cut_cuda.cpp ..
    cp pair_cg_cmm_coul_debye_cuda.cpp ..
    cp pair_cg_cmm_cuda.h ..
    cp pair_cg_cmm_coul_cut_cuda.h ..
    cp pair_cg_cmm_coul_debye_cuda.h ..
  fi

  if (test -e ../pair_cg_cmm_coul_long.cpp) then
    cp pair_cg_cmm_coul_long_cuda.cpp ..
    cp pair_cg_cmm_coul_long_cuda.h ..
  fi

  if (test -e ../pppm.cpp) then
    cp pppm_cuda.cpp ..
    cp fft3d_cuda.cpp ..
    cp fft3d_wrap_cuda.cpp ..
    cp pppm_cuda.h ..
    cp fft3d_cuda.h ..
    cp fft3d_wrap_cuda.h ..
    cp pair_lj_cut_coul_long_cuda.cpp ..
    cp pair_lj_cut_coul_long_cuda.h ..
  fi
  

  if (test -e ../pair_eam.cpp) then
    cp pair_eam_alloy_cuda.cpp ..
    cp pair_eam_cuda.cpp ..
    cp pair_eam_fs_cuda.cpp ..
    cp pair_eam_alloy_cuda.h ..
    cp pair_eam_cuda.h ..
    cp pair_eam_fs_cuda.h ..
  fi

  if (test -e ../pair_gran_hooke.cpp) then
    cp pair_gran_hooke_cuda.cpp ..
    cp pair_gran_hooke_cuda.h ..
  fi

  if (test -e ../pair_lj_charmm_coul_charmm.cpp) then
    cp pair_lj_charmm_coul_charmm_cuda.cpp ..
    cp pair_lj_charmm_coul_charmm_implicit_cuda.cpp ..
    cp pair_lj_charmm_coul_charmm_cuda.h ..
    cp pair_lj_charmm_coul_charmm_implicit_cuda.h ..
    if (test -e ../pair_lj_charmm_coul_long.cpp) then  
      cp pair_lj_charmm_coul_long_cuda.cpp ..
      cp pair_lj_charmm_coul_long_cuda.h ..
    fi
  fi

  if (test -e ../pair_lj_class2.cpp) then
    cp pair_lj_class2_coul_cut_cuda.cpp ..
    cp pair_lj_class2_cuda.cpp ..
    cp pair_lj_class2_coul_cut_cuda.h ..
    cp pair_lj_class2_cuda.h ..
    if (test -e ../pair_lj_class2_coul_long.cpp) then
      cp pair_lj_class2_coul_long_cuda.cpp ..
      cp pair_lj_class2_coul_long_cuda.h ..
    fi
  fi 

  cp atom_vec_atomic_cuda.cpp ..
  cp atom_vec_charge_cuda.cpp ..
  cp comm_cuda.cpp ..
  cp compute_pe_cuda.cpp ..
  cp compute_pressure_cuda.cpp ..
  cp compute_temp_cuda.cpp ..
  cp compute_temp_partial_cuda.cpp ..
  cp domain_cuda.cpp ..
  cp fix_addforce_cuda.cpp ..
  cp fix_aveforce_cuda.cpp ..
  cp fix_enforce2d_cuda.cpp ..
  cp fix_gravity_cuda.cpp ..
  cp fix_nh_cuda.cpp ..
  cp fix_npt_cuda.cpp ..
  cp fix_nve_cuda.cpp ..
  cp fix_nvt_cuda.cpp ..
  cp fix_set_force_cuda.cpp ..
  cp fix_shake_cuda.cpp ..
  cp fix_temp_berendsen_cuda.cpp ..
  cp fix_temp_rescale_cuda.cpp ..
  cp fix_temp_rescale_limit_cuda.cpp ..
  cp fix_viscous_cuda.cpp ..
  cp modify_cuda.cpp ..
  cp neighbor_cuda.cpp ..
  cp neigh_full_cuda.cpp ..
  cp pair_buck_coul_cut_cuda.cpp ..
  cp pair_buck_cuda.cpp ..
  cp pair_lj96_cut_cuda.cpp ..
  cp pair_lj_cut_coul_cut_cuda.cpp ..
  cp pair_lj_cut_coul_debye_cuda.cpp ..
  cp pair_lj_cut_cuda.cpp ..
  cp pair_lj_cut_experimental_cuda.cpp ..
  cp pair_lj_expand_cuda.cpp ..
  cp pair_lj_gromacs_coul_gromacs_cuda.cpp ..
  cp pair_lj_gromacs_cuda.cpp ..
  cp pair_lj_smooth_cuda.cpp ..
  cp pair_morse_cuda.cpp ..
  cp pppm_cuda.cpp ..
  cp verlet_cuda.cpp ..

  cp cuda.cpp ..
  cp cuda_neigh_list.cpp ..

  cp atom_vec_atomic_cuda.h ..
  cp atom_vec_charge_cuda.h ..
  cp comm_cuda.h ..
  cp compute_pe_cuda.h ..
  cp compute_pressure_cuda.h ..
  cp compute_temp_cuda.h ..
  cp compute_temp_partial_cuda.h ..
  cp domain_cuda.h ..
  cp fix_addforce_cuda.h ..
  cp fix_aveforce_cuda.h ..
  cp fix_enforce2d_cuda.h ..
  cp fix_gravity_cuda.h ..
  cp fix_nh_cuda.h ..
  cp fix_npt_cuda.h ..
  cp fix_nve_cuda.h ..
  cp fix_nvt_cuda.h ..
  cp fix_set_force_cuda.h ..
  cp fix_shake_cuda.h ..
  cp fix_temp_berendsen_cuda.h ..
  cp fix_temp_rescale_cuda.h ..
  cp fix_temp_rescale_limit_cuda.h ..
  cp fix_viscous_cuda.h ..
  cp modify_cuda.h ..
  cp neighbor_cuda.h ..
  cp pair_buck_coul_cut_cuda.h ..
  cp pair_buck_cuda.h ..

  cp pair_lj96_cut_cuda.h ..
  cp pair_lj_cut_coul_cut_cuda.h ..
  cp pair_lj_cut_coul_debye_cuda.h ..
  cp pair_lj_cut_cuda.h ..
  cp pair_lj_cut_experimental_cuda.h ..
  cp pair_lj_expand_cuda.h ..
  cp pair_lj_gromacs_coul_gromacs_cuda.h ..
  cp pair_lj_gromacs_cuda.h ..
  cp pair_lj_smooth_cuda.h ..
  cp pair_morse_cuda.h ..
  cp verlet_cuda.h ..

  cp cuda.h ..
  cp cuda_common.h ..
  cp cuda_data.h ..
  cp cuda_modify_flags.h ..
  cp cuda_neigh_list.h ..
  cp cuda_precision.h ..
  cp cuda_shared.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e '/^include.*cuda.*$/d' ../Makefile.package
    sed -i -e 's/[^ \t]*cuda[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*CUDA[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*lrt[^ \t]* //g' ../Makefile.package
  fi

  # force rebuild of dependencies since removing -DLMP_USER_CUDA switch

  touch ../accelerator_cuda.h

  rm ../atom_vec_angle_cuda.cpp
  rm ../atom_vec_atomic_cuda.cpp
  rm ../atom_vec_charge_cuda.cpp
  rm ../atom_vec_full_cuda.cpp
  rm ../comm_cuda.cpp
  rm ../compute_pe_cuda.cpp
  rm ../compute_pressure_cuda.cpp
  rm ../compute_temp_cuda.cpp
  rm ../compute_temp_partial_cuda.cpp
  rm ../domain_cuda.cpp
  rm ../fft3d_cuda.cpp
  rm ../fft3d_wrap_cuda.cpp
  rm ../fix_addforce_cuda.cpp
  rm ../fix_aveforce_cuda.cpp
  rm ../fix_enforce2d_cuda.cpp
  rm ../fix_freeze_cuda.cpp
  rm ../fix_gravity_cuda.cpp
  rm ../fix_nh_cuda.cpp
  rm ../fix_npt_cuda.cpp
  rm ../fix_nve_cuda.cpp
  rm ../fix_nvt_cuda.cpp
  rm ../fix_set_force_cuda.cpp
  rm ../fix_shake_cuda.cpp
  rm ../fix_temp_berendsen_cuda.cpp
  rm ../fix_temp_rescale_cuda.cpp
  rm ../fix_temp_rescale_limit_cuda.cpp
  rm ../fix_viscous_cuda.cpp
  rm ../modify_cuda.cpp
  rm ../neighbor_cuda.cpp
  rm ../neigh_full_cuda.cpp
  rm ../pair_born_coul_long_cuda.cpp
  rm ../pair_buck_coul_cut_cuda.cpp
  rm ../pair_buck_coul_long_cuda.cpp
  rm ../pair_buck_cuda.cpp
  rm ../pair_cg_cmm_coul_cut_cuda.cpp
  rm ../pair_cg_cmm_coul_debye_cuda.cpp
  rm ../pair_cg_cmm_coul_long_cuda.cpp
  rm ../pair_cg_cmm_cuda.cpp
  rm ../pair_eam_alloy_cuda.cpp
  rm ../pair_eam_cuda.cpp
  rm ../pair_eam_fs_cuda.cpp
  rm ../pair_gran_hooke_cuda.cpp
  rm ../pair_lj96_cut_cuda.cpp
  rm ../pair_lj_charmm_coul_charmm_cuda.cpp
  rm ../pair_lj_charmm_coul_charmm_implicit_cuda.cpp
  rm ../pair_lj_charmm_coul_long_cuda.cpp
  rm ../pair_lj_class2_coul_cut_cuda.cpp
  rm ../pair_lj_class2_coul_long_cuda.cpp
  rm ../pair_lj_class2_cuda.cpp
  rm ../pair_lj_cut_coul_cut_cuda.cpp
  rm ../pair_lj_cut_coul_debye_cuda.cpp
  rm ../pair_lj_cut_coul_long_cuda.cpp
  rm ../pair_lj_cut_cuda.cpp
  rm ../pair_lj_cut_experimental_cuda.cpp
  rm ../pair_lj_expand_cuda.cpp
  rm ../pair_lj_gromacs_coul_gromacs_cuda.cpp
  rm ../pair_lj_gromacs_cuda.cpp
  rm ../pair_lj_smooth_cuda.cpp
  rm ../pair_morse_cuda.cpp
  rm ../pppm_cuda.cpp
  rm ../verlet_cuda.cpp

  rm ../cuda.cpp
  rm ../cuda_neigh_list.cpp

  rm ../atom_vec_angle_cuda.h
  rm ../atom_vec_atomic_cuda.h
  rm ../atom_vec_charge_cuda.h
  rm ../atom_vec_full_cuda.h
  rm ../comm_cuda.h
  rm ../compute_pe_cuda.h
  rm ../compute_pressure_cuda.h
  rm ../compute_temp_cuda.h
  rm ../compute_temp_partial_cuda.h
  rm ../domain_cuda.h
  rm ../fft3d_cuda.h
  rm ../fft3d_wrap_cuda.h
  rm ../fix_addforce_cuda.h
  rm ../fix_aveforce_cuda.h
  rm ../fix_enforce2d_cuda.h
  rm ../fix_freeze_cuda.h
  rm ../fix_gravity_cuda.h
  rm ../fix_nh_cuda.h
  rm ../fix_npt_cuda.h
  rm ../fix_nve_cuda.h
  rm ../fix_nvt_cuda.h
  rm ../fix_set_force_cuda.h
  rm ../fix_shake_cuda.h
  rm ../fix_temp_berendsen_cuda.h
  rm ../fix_temp_rescale_cuda.h
  rm ../fix_temp_rescale_limit_cuda.h
  rm ../fix_viscous_cuda.h
  rm ../modify_cuda.h
  rm ../neighbor_cuda.h
  rm ../pair_born_coul_long_cuda.h
  rm ../pair_buck_coul_cut_cuda.h
  rm ../pair_buck_coul_long_cuda.h
  rm ../pair_buck_cuda.h
  rm ../pair_cg_cmm_coul_cut_cuda.h
  rm ../pair_cg_cmm_coul_debye_cuda.h
  rm ../pair_cg_cmm_coul_long_cuda.h
  rm ../pair_cg_cmm_cuda.h
  rm ../pair_eam_alloy_cuda.h
  rm ../pair_eam_cuda.h
  rm ../pair_eam_fs_cuda.h
  rm ../pair_gran_hooke_cuda.h
  rm ../pair_lj96_cut_cuda.h
  rm ../pair_lj_charmm_coul_charmm_cuda.h
  rm ../pair_lj_charmm_coul_charmm_implicit_cuda.h
  rm ../pair_lj_charmm_coul_long_cuda.h
  rm ../pair_lj_class2_coul_cut_cuda.h
  rm ../pair_lj_class2_coul_long_cuda.h
  rm ../pair_lj_class2_cuda.h
  rm ../pair_lj_cut_coul_cut_cuda.h
  rm ../pair_lj_cut_coul_debye_cuda.h
  rm ../pair_lj_cut_coul_long_cuda.h
  rm ../pair_lj_cut_cuda.h
  rm ../pair_lj_cut_experimental_cuda.h
  rm ../pair_lj_expand_cuda.h
  rm ../pair_lj_gromacs_coul_gromacs_cuda.h
  rm ../pair_lj_gromacs_cuda.h
  rm ../pair_lj_smooth_cuda.h
  rm ../pair_morse_cuda.h
  rm ../pppm_cuda.h
  rm ../verlet_cuda.h

  rm ../cuda.h
  rm ../cuda_common.h
  rm ../cuda_data.h
  rm ../cuda_modify_flags.h
  rm ../cuda_neigh_list.h
  rm ../cuda_precision.h
  rm ../cuda_shared.h
fi
