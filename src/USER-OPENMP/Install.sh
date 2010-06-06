# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp dihedral_omp.h ..
  cp dihedral_omp.cpp ..
  
  if (test -e dihedral_harmonic.cpp); then
    cp dihedral_harmonic_omp.h ..
    cp dihedral_harmonic_omp.cpp ..
    cp dihedral_helix_omp.h ..
    cp dihedral_helix_omp.cpp ..
  fi

  cp pair_omp.h ..
  cp pair_omp.cpp ..

  cp pair_buck_omp.h ..
  cp pair_buck_omp.cpp ..
  cp pair_buck_coul_cut_omp.h ..
  cp pair_buck_coul_cut_omp.cpp ..
  cp pair_coul_cut_omp.h ..
  cp pair_coul_cut_omp.cpp ..
  cp pair_coul_debye_omp.h ..
  cp pair_coul_debye_omp.cpp ..
  cp pair_dpd_omp.h ..
  cp pair_dpd_omp.cpp ..
  cp pair_dpd_tstat_omp.h ..
  cp pair_dpd_tstat_omp.cpp ..
  cp pair_lj_charmm_coul_charmm_implicit_omp.h ..
  cp pair_lj_charmm_coul_charmm_implicit_omp.cpp ..
  cp pair_lj_charmm_coul_charmm_omp.h ..
  cp pair_lj_charmm_coul_charmm_omp.cpp ..
  cp pair_lj_cut_coul_cut_omp.h ..
  cp pair_lj_cut_coul_cut_omp.cpp ..
  cp pair_lj_cut_coul_debye_omp.h ..
  cp pair_lj_cut_coul_debye_omp.cpp ..
  cp pair_lj_cut_omp.h ..
  cp pair_lj_cut_omp.cpp ..
  cp pair_lj_expand_omp.h ..
  cp pair_lj_expand_omp.cpp ..
  cp pair_lj_gromacs_coul_gromacs_omp.h ..
  cp pair_lj_gromacs_coul_gromacs_omp.cpp ..
  cp pair_lj96_cut_omp.h ..
  cp pair_lj96_cut_omp.cpp ..
  cp pair_lj_gromacs_omp.h ..
  cp pair_lj_gromacs_omp.cpp ..
  cp pair_lj_smooth_omp.h ..
  cp pair_lj_smooth_omp.cpp ..
  cp pair_morse_omp.h ..
  cp pair_morse_omp.cpp ..
  cp pair_soft_omp.h ..
  cp pair_soft_omp.cpp ..
  cp pair_table_omp.h ..
  cp pair_table_omp.cpp ..
  cp pair_yukawa_omp.h ..
  cp pair_yukawa_omp.cpp ..

  if (test -e ../pair_cg_cmm.cpp) then
    cp pair_cg_cmm_omp.h ..
    cp pair_cg_cmm_omp.cpp ..
  fi

  if (test -e pair_colloid.cpp); then
    cp pair_colloid_omp.h ..
    cp pair_colloid_omp.cpp ..
    cp pair_lubricate_omp.h ..
    cp pair_lubricate_omp.cpp ..
    cp pair_yukawa_colloid_omp.h ..
    cp pair_yukawa_colloid_omp.cpp ..
  fi

  if (test -e pair_dipole_cut.cpp); then
    cp pair_dipole_cut_omp.h ..
    cp pair_dipole_cut_omp.cpp ..
  fi

  if (test -e pair_gayberne.cpp); then
    cp pair_gayberne_omp.h ..
    cp pair_gayberne_omp.cpp ..
    cp pair_resquared_omp.h ..
    cp pair_resquared_omp.cpp ..
  fi

  if (test -e ../pair_airebo.cpp) then
    cp pair_airebo_omp.h ..
    cp pair_airebo_omp.cpp ..
    cp pair_eam_omp.h ..
    cp pair_eam_omp.cpp ..
    cp pair_eam_alloy_omp.h ..
    cp pair_eam_alloy_omp.cpp ..
    cp pair_eam_fs_omp.h ..
    cp pair_eam_fs_omp.cpp ..
    cp pair_sw_omp.h ..
    cp pair_sw_omp.cpp ..
    cp pair_tersoff_omp.h ..
    cp pair_tersoff_omp.cpp ..
    cp pair_tersoff_zbl_omp.h ..
    cp pair_tersoff_zbl_omp.cpp ..
  fi

  if (test -e ../pair_gauss.cpp) then
    cp pair_coul_diel_omp.h ..
    cp pair_coul_diel_omp.cpp ..
    cp pair_gauss_cut_omp.h ..
    cp pair_gauss_cut_omp.cpp ..
  fi

  if (test -e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_born_coul_long_omp.h ..
    cp pair_born_coul_long_omp.cpp ..
    cp pair_buck_coul_long_omp.h ..
    cp pair_buck_coul_long_omp.cpp ..
    cp pair_coul_long_omp.h ..
    cp pair_coul_long_omp.cpp ..
    cp pair_lj_charmm_coul_long_omp.h ..
    cp pair_lj_charmm_coul_long_omp.cpp ..
    cp pair_lj_cut_coul_long_omp.h ..
    cp pair_lj_cut_coul_long_omp.cpp ..
    cp ewald_omp.h ..
    cp ewald_omp.cpp ..
  fi

elif (test $1 = 0) then

  rm ../dihedral_omp.h
  rm ../dihedral_omp.cpp
  rm -f ../dihedral_harmonic_omp.h
  rm -f ../dihedral_harmonic_omp.cpp
  rm -f ../dihedral_helix_omp.h
  rm -f ../dihedral_helix_omp.cpp

  rm ../pair_omp.h
  rm ../pair_omp.cpp

  rm ../pair_buck_omp.h
  rm ../pair_buck_omp.cpp
  rm ../pair_buck_coul_cut_omp.h
  rm ../pair_buck_coul_cut_omp.cpp
  rm ../pair_coul_cut_omp.h
  rm ../pair_coul_cut_omp.cpp
  rm ../pair_coul_debye_omp.h
  rm ../pair_coul_debye_omp.cpp
  rm ../pair_dpd_omp.h
  rm ../pair_dpd_omp.cpp
  rm ../pair_dpd_tstat_omp.h
  rm ../pair_dpd_tstat_omp.cpp
  rm ../pair_lj_charmm_coul_charmm_implicit_omp.h
  rm ../pair_lj_charmm_coul_charmm_implicit_omp.cpp
  rm ../pair_lj_charmm_coul_charmm_omp.h
  rm ../pair_lj_charmm_coul_charmm_omp.cpp
  rm ../pair_lj_cut_coul_cut_omp.h
  rm ../pair_lj_cut_coul_cut_omp.cpp
  rm ../pair_lj_cut_coul_debye_omp.h
  rm ../pair_lj_cut_coul_debye_omp.cpp
  rm ../pair_lj_cut_omp.h
  rm ../pair_lj_cut_omp.cpp
  rm ../pair_lj_expand_omp.h
  rm ../pair_lj_expand_omp.cpp
  rm ../pair_lj_gromacs_coul_gromacs_omp.h
  rm ../pair_lj_gromacs_coul_gromacs_omp.cpp
  rm ../pair_lj96_cut_omp.h
  rm ../pair_lj96_cut_omp.cpp
  rm ../pair_lj_gromacs_omp.h
  rm ../pair_lj_gromacs_omp.cpp
  rm ../pair_lj_smooth_omp.h
  rm ../pair_lj_smooth_omp.cpp
  rm ../pair_morse_omp.h
  rm ../pair_morse_omp.cpp
  rm ../pair_soft_omp.h
  rm ../pair_soft_omp.cpp
  rm ../pair_table_omp.h
  rm ../pair_table_omp.cpp
  rm ../pair_yukawa_omp.h
  rm ../pair_yukawa_omp.cpp

  rm -f ../pair_cg_cmm_omp.h
  rm -f ../pair_cg_cmm_omp.cpp

  rm -f ../pair_colloid_omp.h
  rm -f ../pair_colloid_omp.cpp
  rm -f ../pair_lubricate_omp.h
  rm -f ../pair_lubricate_omp.cpp
  rm -f ../pair_yukawa_colloid_omp.h
  rm -f ../pair_yukawa_colloid_omp.cpp

  rm -f ../pair_dipole_cut_omp.h
  rm -f ../pair_dipole_cut_omp.cpp

  rm -f ../pair_gayberne_omp.h
  rm -f ../pair_gayberne_omp.cpp
  rm -f ../pair_resquared_omp.h
  rm -f ../pair_resquared_omp.cpp

  rm -f ../pair_airebo_omp.h
  rm -f ../pair_airebo_omp.cpp
  rm -f ../pair_eam_omp.h
  rm -f ../pair_eam_omp.cpp
  rm -f ../pair_eam_alloy_omp.h
  rm -f ../pair_eam_alloy_omp.cpp
  rm -f ../pair_eam_fs_omp.h
  rm -f ../pair_eam_fs_omp.cpp
  rm -f ../pair_sw_omp.h
  rm -f ../pair_sw_omp.cpp
  rm -f ../pair_tersoff_omp.h
  rm -f ../pair_tersoff_omp.cpp
  rm -f ../pair_tersoff_zbl_omp.h
  rm -f ../pair_tersoff_zbl_omp.cpp

  rm -f ../pair_coul_diel_omp.h
  rm -f ../pair_coul_diel_omp.cpp
  rm -f ../pair_gauss_cut_omp.h
  rm -f ../pair_gauss_cut_omp.cpp

  rm -f ../pair_born_coul_long_omp.h
  rm -f ../pair_born_coul_long_omp.cpp
  rm -f ../pair_buck_coul_long_omp.h
  rm -f ../pair_buck_coul_long_omp.cpp
  rm -f ../pair_coul_long_omp.h
  rm -f ../pair_coul_long_omp.cpp
  rm -f ../pair_lj_charmm_coul_long_omp.h
  rm -f ../pair_lj_charmm_coul_long_omp.cpp
  rm -f ../pair_lj_cut_coul_long_omp.h
  rm -f ../pair_lj_cut_coul_long_omp.cpp
  rm -f ../ewald_omp.h
  rm -f ../ewald_omp.cpp
fi
