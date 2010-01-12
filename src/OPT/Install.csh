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

  if (-e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_lj_charmm_coul_long_opt.cpp ..
    cp pair_lj_charmm_coul_long_opt.h ..
  endif

  cp pair_lj_cut_opt.cpp ..
  cp pair_lj_cut_opt.h ..

  cp pair_morse_opt.cpp ..
  cp pair_morse_opt.h ..

else if ($1 == 0) then

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
