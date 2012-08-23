# Update package files in LAMMPS
# cp package file to src if doesn't exist or is different
# do not copy gayberne files if non-GPU version does not exist

for file in *.cpp *.h; do
  if (test $file = pair_gayberne_gpu.cpp -a ! -e ../pair_gayberne.cpp) then
    continue
  fi
  if (test $file = pair_gayberne_gpu.h -a ! -e ../pair_gayberne.cpp) then
    continue
  fi
  if (test $file = pair_resquared_gpu.cpp -a ! -e ../pair_resquared.cpp) then
    continue
  fi
  if (test $file = pair_resquared_gpu.h -a ! -e ../pair_resquared.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_gpu.cpp -a ! -e ../pair_lj_cut_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_gpu.h -a ! -e ../pair_lj_cut_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_coul_long_gpu.cpp -a ! -e ../pair_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_coul_long_gpu.h -a ! -e ../pair_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_sdk_gpu.cpp -a ! -e ../pair_lj_sdk.cpp) then
    continue
  fi
  if (test $file = pair_lj_sdk_gpu.h -a ! -e ../pair_lj_sdk.cpp) then
    continue
  fi
  if (test $file = pair_lj_sdk_coul_long_gpu.cpp -a ! -e ../pair_lj_sdk_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_sdk_coul_long_gpu.h -a ! -e ../pair_lj_sdk_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_dipole_sf_gpu.cpp -a ! -e ../pair_dipole_sf.cpp) then
    continue
  fi
  if (test $file = pair_dipole_sf_gpu.h -a ! -e ../pair_dipole_sf.cpp) then
    continue
  fi
  if (test $file = pair_dipole_cut_gpu.cpp -a ! -e ../pair_dipole_cut.cpp) then
    continue
  fi
  if (test $file = pair_dipole_cut_gpu.h -a ! -e ../pair_dipole_cut.cpp) then
    continue
  fi
  if (test $file = pair_yukawa_colloid_gpu.cpp -a ! -e ../pair_yukawa_colloid.cpp) then
    continue
  fi
  if (test $file = pair_yukawa_colloid_gpu.h -a ! -e ../pair_yukawa_colloid.cpp) then
    continue
  fi
  if (test $file = pair_colloid_gpu.cpp -a ! -e ../pair_colloid.cpp) then
    continue
  fi
  if (test $file = pair_colloid_gpu.h -a ! -e ../pair_colloid.cpp) then
    continue
  fi
  if (test $file = pair_buck_coul_long_gpu.cpp -a ! -e ../pair_buck_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_buck_coul_long_gpu.h -a ! -e ../pair_buck_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_born_coul_long_gpu.cpp -a ! -e ../pair_born_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_born_coul_long_gpu.h -a ! -e ../pair_born_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_eam_gpu.cpp -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_gpu.h -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_alloy_gpu.cpp -a ! -e ../pair_eam_alloy.cpp) then
    continue
  fi
  if (test $file = pair_eam_alloy_gpu.h -a ! -e ../pair_eam_alloy.cpp) then
    continue
  fi
  if (test $file = pair_eam_fs_gpu.cpp -a ! -e ../pair_eam_fs.cpp) then
    continue
  fi
  if (test $file = pair_eam_fs_gpu.h -a ! -e ../pair_eam_fs.cpp) then
    continue
  fi
  if (test $file = pair_lj_class2_gpu.cpp -a ! -e ../pair_lj_class2.cpp) then
    continue
  fi
  if (test $file = pair_lj_class2_coul_long_gpu.cpp -a ! -e ../pair_lj_class2_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_class2_coul_long_gpu.h -a ! -e ../pair_lj_class2_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_charmm_coul_long_gpu.cpp -a ! -e ../pair_lj_charmm_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_charmm_coul_long_gpu.h -a ! -e ../pair_lj_charmm_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_dsf_gpu.cpp -a ! -e ../pair_lj_cut_coul_dsf.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_dsf_gpu.h -a ! -e ../pair_lj_cut_coul_dsf.cpp) then
    continue
  fi
  if (test $file = pppm_gpu.cpp -a ! -e ../pppm.cpp) then
    continue
  fi
  if (test $file = pppm_gpu.h -a ! -e ../pppm.cpp) then
    continue
  fi
  
  if (test ! -e ../$file) then
    echo "  creating src/$file"
    cp $file ..
  elif (test "`diff --brief $file ../$file`" != "") then
    echo "  updating src/$file"
    cp $file ..
  fi
done
