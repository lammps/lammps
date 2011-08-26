# Update package files in LAMMPS
# cp package file to src if doesn't exist or is different
# do not copy eam and charmm files if non-OPT versions do not exist

for file in *.cpp *.h; do
  if (test $file = pair_eam_opt.cpp -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_opt.h -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_alloy_opt.cpp -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_alloy_opt.h -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_fs_opt.cpp -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_eam_fs_opt.h -a ! -e ../pair_eam.cpp) then
    continue
  fi
  if (test $file = pair_lj_charmm_coul_long_opt.cpp -a ! -e ../pair_lj_charmm_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_charmm_coul_long_opt.h -a ! -e ../pair_lj_charmm_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_opt.cpp -a ! -e ../pair_lj_cut_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_opt.h -a ! -e ../pair_lj_cut_coul_long.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_tip4p_opt.cpp -a ! -e ../pair_lj_cut_coul_long_tip4p.cpp) then
    continue
  fi
  if (test $file = pair_lj_cut_coul_long_tip4p_opt.h -a ! -e ../pair_lj_cut_coul_long_tip4p.cpp) then
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
