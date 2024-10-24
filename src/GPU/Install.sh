# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

cat <<EOF
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING

  Support for building the GPU package with the legacy build system using GNU
  make will be removed in Summer 2025.  Please switch to using CMake to build
  LAMMPS as soon as possible and report any problems to developers@lammps.org
  or post a bug report issue at https://github.com/lammps/lammps/issues

WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
EOF

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

# list of files with optional dependcies

action amoeba_convolution_gpu.cpp amoeba_convolution.cpp
action amoeba_convolution_gpu.h amoeba_convolution.cpp
action fix_gpu.cpp
action fix_gpu.h
action fix_nve_gpu.h
action fix_nve_gpu.cpp
action fix_nh_gpu.h
action fix_nh_gpu.cpp
action fix_nvt_gpu.h
action fix_nvt_gpu.cpp
action fix_npt_gpu.h
action fix_npt_gpu.cpp
action fix_nve_asphere_gpu.h fix_nve_asphere.h
action fix_nve_asphere_gpu.cpp fix_nve_asphere.cpp
action gpu_extra.h
action pair_amoeba_gpu.cpp pair_amoeba.cpp
action pair_amoeba_gpu.h pair_amoeba.h
action pair_beck_gpu.cpp pair_beck.cpp
action pair_beck_gpu.h pair_beck.h
action pair_born_coul_long_gpu.cpp pair_born_coul_long.cpp
action pair_born_coul_long_gpu.h pair_born_coul_long.cpp
action pair_born_coul_long_cs_gpu.cpp pair_born_coul_long_cs.cpp
action pair_born_coul_long_cs_gpu.h pair_born_coul_long_cs.cpp
action pair_born_coul_wolf_gpu.cpp pair_born_coul_wolf.cpp
action pair_born_coul_wolf_gpu.h pair_born_coul_wolf.h
action pair_born_coul_wolf_cs_gpu.cpp pair_born_coul_wolf_cs.cpp
action pair_born_coul_wolf_cs_gpu.h pair_born_coul_wolf_cs.cpp
action pair_born_gpu.cpp
action pair_born_gpu.h
action pair_buck_coul_cut_gpu.cpp pair_buck_coul_cut.cpp
action pair_buck_coul_cut_gpu.h pair_buck_coul_cut.cpp
action pair_buck_coul_long_gpu.cpp pair_buck_coul_long.cpp
action pair_buck_coul_long_gpu.h pair_buck_coul_long.cpp
action pair_buck_gpu.cpp pair_buck.cpp
action pair_buck_gpu.h pair_buck.cpp
action pair_colloid_gpu.cpp pair_colloid.cpp
action pair_colloid_gpu.h pair_colloid.cpp
action pair_coul_cut_gpu.cpp
action pair_coul_cut_gpu.h
action pair_coul_debye_gpu.cpp
action pair_coul_debye_gpu.h
action pair_coul_dsf_gpu.cpp
action pair_coul_dsf_gpu.h
action pair_coul_long_gpu.cpp pair_coul_long.cpp
action pair_coul_long_gpu.h pair_coul_long.cpp
action pair_coul_long_cs_gpu.cpp pair_coul_long_cs.cpp
action pair_coul_long_cs_gpu.h pair_coul_long_cs.cpp
action pair_dpd_gpu.cpp pair_dpd.cpp
action pair_dpd_gpu.h pair_dpd.h
action pair_dpd_tstat_gpu.cpp pair_dpd_tstat.cpp
action pair_dpd_tstat_gpu.h pair_dpd_tstat.h
action pair_lj_cut_dipole_cut_gpu.cpp pair_lj_cut_dipole_cut.cpp
action pair_lj_cut_dipole_cut_gpu.h pair_lj_cut_dipole_cut.cpp
action pair_lj_sf_dipole_sf_gpu.cpp pair_lj_sf_dipole_sf.cpp
action pair_lj_sf_dipole_sf_gpu.h pair_lj_sf_dipole_sf.cpp
action pair_eam_alloy_gpu.cpp pair_eam.cpp
action pair_eam_alloy_gpu.h pair_eam.cpp
action pair_eam_fs_gpu.cpp pair_eam.cpp
action pair_eam_fs_gpu.h pair_eam.cpp
action pair_eam_gpu.cpp pair_eam.cpp
action pair_eam_gpu.h pair_eam.cpp
action pair_gauss_gpu.cpp pair_gauss.cpp
action pair_gauss_gpu.h pair_gauss.h
action pair_gayberne_gpu.cpp pair_gayberne.cpp
action pair_gayberne_gpu.h pair_gayberne.cpp
action pair_hippo_gpu.cpp pair_hippo.cpp
action pair_hippo_gpu.h pair_hippo.cpp
action pair_lj96_cut_gpu.cpp pair_lj96_cut.cpp
action pair_lj96_cut_gpu.h pair_lj96_cut.h
action pair_lj_charmm_coul_long_gpu.cpp pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_long_gpu.h pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_charmm_gpu.cpp pair_lj_charmm_coul_charmm.cpp
action pair_lj_charmm_coul_charmm_gpu.h pair_lj_charmm_coul_charmm.cpp
action pair_lj_class2_coul_long_gpu.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_gpu.h pair_lj_class2_coul_long.cpp
action pair_lj_class2_gpu.cpp pair_lj_class2.cpp
action pair_lj_class2_gpu.h pair_lj_class2.cpp
action pair_lj_cubic_gpu.cpp pair_lj_cubic.cpp
action pair_lj_cubic_gpu.h pair_lj_cubic.h
action pair_lj_cut_coul_cut_gpu.cpp
action pair_lj_cut_coul_cut_gpu.h
action pair_lj_cut_coul_debye_gpu.cpp pair_lj_cut_coul_debye.cpp
action pair_lj_cut_coul_debye_gpu.h pair_lj_cut_coul_debye.h
action pair_lj_cut_coul_dsf_gpu.cpp pair_lj_cut_coul_dsf.cpp
action pair_lj_cut_coul_dsf_gpu.h pair_lj_cut_coul_dsf.h
action pair_lj_cut_coul_long_gpu.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_gpu.h pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_msm_gpu.cpp pair_lj_cut_coul_msm.cpp
action pair_lj_cut_coul_msm_gpu.h pair_lj_cut_coul_msm.h
action pair_lj_cut_gpu.cpp
action pair_lj_cut_gpu.h
action pair_lj_cut_dipole_long_gpu.cpp pair_lj_cut_dipole_long.cpp
action pair_lj_cut_dipole_long_gpu.h pair_lj_cut_dipole_long.cpp
action pair_lj_cut_tip4p_long_gpu.h pair_lj_cut_tip4p_long.cpp
action pair_lj_cut_tip4p_long_gpu.cpp pair_lj_cut_tip4p_long.cpp
action pair_lj_smooth_gpu.cpp pair_lj_smooth.cpp
action pair_lj_smooth_gpu.h pair_lj_smooth.cpp
action pair_lj_expand_gpu.cpp
action pair_lj_expand_gpu.h
action pair_lj_expand_coul_long_gpu.cpp pair_lj_expand_coul_long.cpp
action pair_lj_expand_coul_long_gpu.h pair_lj_expand_coul_long.cpp
action pair_lj_gromacs_gpu.cpp pair_lj_gromacs.cpp
action pair_lj_gromacs_gpu.h pair_lj_gromacs.h
action pair_lj_spica_coul_long_gpu.cpp pair_lj_spica_coul_long.cpp
action pair_lj_spica_coul_long_gpu.h pair_lj_spica_coul_long.cpp
action pair_lj_spica_gpu.cpp pair_lj_spica.cpp
action pair_lj_spica_gpu.h pair_lj_spica.cpp
action pair_mie_cut_gpu.cpp pair_mie_cut.cpp
action pair_mie_cut_gpu.h pair_mie_cut.h
action pair_morse_gpu.cpp
action pair_morse_gpu.h
action pair_resquared_gpu.cpp pair_resquared.cpp
action pair_resquared_gpu.h pair_resquared.cpp
action pair_soft_gpu.cpp
action pair_soft_gpu.h
action pair_sw_gpu.cpp pair_sw.cpp
action pair_sw_gpu.h pair_sw.h
action pair_vashishta_gpu.cpp pair_vashishta.cpp
action pair_vashishta_gpu.h pair_vashishta.h
action pair_table_gpu.cpp pair_table.cpp
action pair_table_gpu.h pair_table.cpp
action pair_tersoff_gpu.cpp pair_tersoff.cpp
action pair_tersoff_gpu.h pair_tersoff.cpp
action pair_tersoff_mod_gpu.cpp pair_tersoff_mod.cpp
action pair_tersoff_mod_gpu.h pair_tersoff_mod.cpp
action pair_tersoff_zbl_gpu.cpp pair_tersoff_zbl.cpp
action pair_tersoff_zbl_gpu.h pair_tersoff_zbl.cpp
action pair_yukawa_colloid_gpu.cpp pair_yukawa_colloid.cpp
action pair_yukawa_colloid_gpu.h pair_yukawa_colloid.cpp
action pair_yukawa_gpu.cpp pair_yukawa.cpp
action pair_yukawa_gpu.h pair_yukawa.cpp
action pair_zbl_gpu.cpp
action pair_zbl_gpu.h
action pppm_gpu.cpp pppm.cpp
action pppm_gpu.h pppm.cpp
action pair_ufm_gpu.cpp pair_ufm.cpp
action pair_ufm_gpu.h pair_ufm.h

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*GPU[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_GPU |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(gpu_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*gpu.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/gpu\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*GPU[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*gpu.*$/d' ../Makefile.package.settings
  fi

fi
