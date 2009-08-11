# Install/unInstall package classes in LAMMPS
# do not copy eam and charmm files if non-OPT versions do not exist

if ($1 == 1) then

  if (-e ../pair_eam.cpp) then
    cp pair_eam_opt.cpp ..
    cp pair_eam_alloy_opt.cpp ..
    cp pair_eam_fs_opt.cpp ..
    cp pair_eam_opt.h ..
    cp pair_eam_alloy_opt.h ..
    cp pair_eam_fs_opt.h ..
  endif

  cp style_opt.h tmp.h
  if (! -e ../pair_eam.cpp) then
    grep -v eam tmp.h > tmp1.h
    mv tmp1.h tmp.h
  endif
  if (! -e ../pair_lj_charmm_coul_long.cpp) then
    grep -v charmm tmp.h > tmp1.h
    mv tmp1.h tmp.h
  endif
  mv tmp.h ../style_opt.h

  if (-e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_lj_charmm_coul_long_opt.cpp ..
    cp pair_lj_charmm_coul_long_opt.h ..
  endif

  cp pair_lj_cut_opt.cpp ..
  cp pair_lj_cut_opt.h ..

  cp pair_morse_opt.cpp ..
  cp pair_morse_opt.h ..

else if ($1 == 0) then

  rm ../style_opt.h
  touch ../style_opt.h

  rm ../pair_eam_opt.cpp
  rm ../pair_eam_alloy_opt.cpp
  rm ../pair_eam_fs_opt.cpp
  rm ../pair_lj_charmm_coul_long_opt.cpp
  rm ../pair_lj_cut_opt.cpp
  rm ../pair_morse_opt.cpp

  rm ../pair_eam_opt.h
  rm ../pair_eam_alloy_opt.h
  rm ../pair_eam_fs_opt.h
  rm ../pair_lj_charmm_coul_long_opt.h
  rm ../pair_lj_cut_opt.h
  rm ../pair_morse_opt.h

endif
