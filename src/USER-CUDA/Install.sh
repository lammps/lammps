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

# force rebuild of files with LMP_USER_CUDA switch

touch ../accelerator_cuda.h

# list of files with optional dependencies

action atom_vec_angle_cuda.cpp atom_vec_angle.cpp
action atom_vec_angle_cuda.h atom_vec_angle.cpp
action atom_vec_atomic_cuda.cpp
action atom_vec_atomic_cuda.h
action atom_vec_charge_cuda.cpp
action atom_vec_charge_cuda.h
action atom_vec_full_cuda.cpp atom_vec_full.cpp
action atom_vec_full_cuda.h atom_vec_full.cpp
action comm_cuda.cpp
action comm_cuda.h
action compute_pe_cuda.cpp
action compute_pe_cuda.h
action compute_pressure_cuda.cpp
action compute_pressure_cuda.h
action compute_temp_cuda.cpp
action compute_temp_cuda.h
action compute_temp_partial_cuda.cpp
action compute_temp_partial_cuda.h
action cuda.cpp
action cuda.h
action cuda_data.h
action cuda_modify_flags.h
action cuda_neigh_list.cpp
action cuda_neigh_list.h
action domain_cuda.cpp
action domain_cuda.h
action fft3d_cuda.cpp pppm.cpp
action fft3d_cuda.h pppm.cpp
action fft3d_wrap_cuda.cpp pppm.cpp
action fft3d_wrap_cuda.h pppm.cpp
action fix_addforce_cuda.cpp
action fix_addforce_cuda.h
action fix_aveforce_cuda.cpp
action fix_aveforce_cuda.h
action fix_enforce2d_cuda.cpp
action fix_enforce2d_cuda.h
action fix_freeze_cuda.cpp fix_freeze.cpp
action fix_freeze_cuda.h fix_freeze.cpp
action fix_gravity_cuda.cpp
action fix_gravity_cuda.h
action fix_nh_cuda.cpp
action fix_nh_cuda.h
action fix_npt_cuda.cpp
action fix_npt_cuda.h
action fix_nve_cuda.cpp
action fix_nve_cuda.h
action fix_nvt_cuda.cpp
action fix_nvt_cuda.h
action fix_set_force_cuda.cpp
action fix_set_force_cuda.h
action fix_shake_cuda.cpp
action fix_shake_cuda.h
action fix_temp_berendsen_cuda.cpp
action fix_temp_berendsen_cuda.h
action fix_temp_rescale_cuda.cpp
action fix_temp_rescale_cuda.h
action fix_temp_rescale_limit_cuda.cpp
action fix_temp_rescale_limit_cuda.h
action fix_viscous_cuda.cpp
action fix_viscous_cuda.h
action modify_cuda.cpp
action modify_cuda.h
action neigh_full_cuda.cpp
action neighbor_cuda.cpp
action neighbor_cuda.h
action pair_born_coul_long_cuda.cpp pair_born_coul_long.cpp
action pair_born_coul_long_cuda.h pair_born_coul_long.cpp
action pair_buck_coul_cut_cuda.cpp
action pair_buck_coul_cut_cuda.h
action pair_buck_coul_long_cuda.cpp pair_buck_coul_long.cpp
action pair_buck_coul_long_cuda.h pair_buck_coul_long.cpp
action pair_buck_cuda.cpp
action pair_buck_cuda.h
action pair_eam_alloy_cuda.cpp pair_eam_alloy.cpp
action pair_eam_alloy_cuda.h pair_eam_alloy.cpp
action pair_eam_cuda.cpp pair_eam.cpp
action pair_eam_cuda.h pair_eam.cpp
action pair_eam_fs_cuda.cpp pair_eam_fs.cpp
action pair_eam_fs_cuda.h pair_eam_fs.cpp
action pair_gran_hooke_cuda.cpp pair_gran_hooke.cpp
action pair_gran_hooke_cuda.h pair_gran_hooke.cpp
action pair_lj96_cut_cuda.cpp
action pair_lj96_cut_cuda.h
action pair_lj_charmm_coul_charmm_cuda.cpp pair_lj_charmm_coul_charmm.cpp
action pair_lj_charmm_coul_charmm_cuda.h pair_lj_charmm_coul_charmm.cpp
action pair_lj_charmm_coul_charmm_implicit_cuda.cpp pair_lj_charmm_coul_charmm_implicit.cpp
action pair_lj_charmm_coul_charmm_implicit_cuda.h pair_lj_charmm_coul_charmm_implicit.cpp
action pair_lj_charmm_coul_long_cuda.cpp pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_long_cuda.h pair_lj_charmm_coul_long.cpp
action pair_lj_class2_coul_cut_cuda.cpp pair_lj_class2_coul_cut.cpp
action pair_lj_class2_coul_cut_cuda.h pair_lj_class2_coul_cut.cpp
action pair_lj_class2_coul_long_cuda.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_cuda.h pair_lj_class2_coul_long.cpp
action pair_lj_class2_cuda.cpp pair_lj_class2.cpp
action pair_lj_class2_cuda.h pair_lj_class2.cpp
action pair_lj_cut_coul_cut_cuda.cpp
action pair_lj_cut_coul_cut_cuda.h
action pair_lj_cut_coul_debye_cuda.cpp
action pair_lj_cut_coul_debye_cuda.h
action pair_lj_cut_coul_long_cuda.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_cuda.h pair_lj_cut_coul_long.cpp
action pair_lj_cut_cuda.cpp
action pair_lj_cut_cuda.h
action pair_lj_cut_experimental_cuda.cpp
action pair_lj_cut_experimental_cuda.h
action pair_lj_expand_cuda.cpp
action pair_lj_expand_cuda.h
action pair_lj_gromacs_coul_gromacs_cuda.cpp
action pair_lj_gromacs_coul_gromacs_cuda.h
action pair_lj_gromacs_cuda.cpp
action pair_lj_gromacs_cuda.h
action pair_lj_sdk_coul_long_cuda.cpp pair_lj_sdk_coul_long.cpp
action pair_lj_sdk_coul_long_cuda.h pair_lj_sdk_coul_long.cpp
action pair_lj_sdk_cuda.cpp pair_lj_sdk.cpp
action pair_lj_sdk_cuda.h pair_lj_sdk.cpp
action pair_lj_smooth_cuda.cpp
action pair_lj_smooth_cuda.h
action pair_morse_cuda.cpp
action pair_morse_cuda.h
action pair_sw_cuda.cpp pair_sw.cpp
action pair_sw_cuda.h pair_sw.cpp
action pair_tersoff_cuda.cpp pair_tersoff.cpp
action pair_tersoff_cuda.h pair_tersoff.cpp
action pair_tersoff_zbl_cuda.cpp pair_tersoff_zbl.cpp
action pair_tersoff_zbl_cuda.h pair_tersoff_zbl.cpp
action pppm_cuda.cpp pppm.cpp
action pppm_cuda.h pppm.cpp
action verlet_cuda.cpp
action verlet_cuda.h

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*cuda[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*CUDA[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I..\/..\/lib\/cuda -DLMP_USER_CUDA |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/cuda |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-llammpscuda |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(user-cuda_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-cuda_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-cuda_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*cuda.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/cuda\/Makefile.lammps
' ../Makefile.package.settings

  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*cuda[^ \t]* //g' ../Makefile.package
    sed -i -e 's/[^ \t]*CUDA[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*cuda.*$/d' ../Makefile.package.settings
  fi

fi
