# Install/unInstall package files in LAMMPS
# do not install child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../pair_eam.cpp) then
    cp pair_eam_opt.cpp ..
    cp pair_eam_alloy_opt.cpp ..
    cp pair_eam_fs_opt.cpp ..
    cp pair_eam_opt.h ..
    cp pair_eam_alloy_opt.h ..
    cp pair_eam_fs_opt.h ..
  fi

  if (test -e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_lj_charmm_coul_long_opt.cpp ..
    cp pair_lj_charmm_coul_long_opt.h ..
    cp pair_lj_cut_coul_long_opt.cpp ..
    cp pair_lj_cut_coul_long_opt.h ..
    cp pair_lj_cut_coul_long_tip4p_opt.cpp ..
    cp pair_lj_cut_coul_long_tip4p_opt.h ..
  fi

  cp pair_lj_cut_opt.cpp ..
  cp pair_lj_cut_opt.h ..

  cp pair_morse_opt.cpp ..
  cp pair_morse_opt.h ..

elif (test $1 = 0) then

  rm -f ../pair_eam_opt.cpp
  rm -f ../pair_eam_alloy_opt.cpp
  rm -f ../pair_eam_fs_opt.cpp
  rm -f ../pair_lj_charmm_coul_long_opt.cpp
  rm -f ../pair_lj_cut_coul_long_opt.cpp
  rm -f ../pair_lj_cut_coul_long_tip4p_opt.cpp
  rm -f ../pair_lj_cut_opt.cpp
  rm -f ../pair_morse_opt.cpp

  rm -f ../pair_eam_opt.h
  rm -f ../pair_eam_alloy_opt.h
  rm -f ../pair_eam_fs_opt.h
  rm -f ../pair_lj_charmm_coul_long_opt.h
  rm -f ../pair_lj_cut_coul_long_opt.h
  rm -f ../pair_lj_cut_coul_long_tip4p_opt.h
  rm -f ../pair_lj_cut_opt.h
  rm -f ../pair_morse_opt.h

fi
