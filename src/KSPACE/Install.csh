# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_kspace.h ..

  cp ewald.cpp ..
  cp pppm.cpp ..
  cp pppm_tip4p.cpp ..
  cp pair_born_coul_long.cpp ..
  cp pair_buck_coul_long.cpp ..
  cp pair_coul_long.cpp ..
  cp pair_lj_cut_coul_long.cpp ..
  cp pair_lj_cut_coul_long_tip4p.cpp ..
  cp pair_lj_charmm_coul_long.cpp ..
  cp fft3d.cpp ..
  cp fft3d_wrap.cpp ..
  cp remap.cpp ..
  cp remap_wrap.cpp ..

  cp ewald.h ..
  cp pppm.h ..
  cp pppm_tip4p.h ..
  cp pair_born_coul_long.h ..
  cp pair_buck_coul_long.h ..
  cp pair_coul_long.h ..
  cp pair_lj_cut_coul_long.h ..
  cp pair_lj_cut_coul_long_tip4p.h ..
  cp pair_lj_charmm_coul_long.h ..
  cp fft3d.h ..
  cp fft3d_wrap.h ..
  cp remap.h ..
  cp remap_wrap.h ..

else if ($1 == 0) then

  rm ../style_kspace.h
  touch ../style_kspace.h

  rm ../ewald.cpp
  rm ../pppm.cpp
  rm ../pppm_tip4p.cpp
  rm ../pair_born_coul_long.cpp
  rm ../pair_buck_coul_long.cpp
  rm ../pair_coul_long.cpp
  rm ../pair_lj_cut_coul_long.cpp
  rm ../pair_lj_cut_coul_long_tip4p.cpp
  rm ../pair_lj_charmm_coul_long.cpp
  rm ../fft3d.cpp
  rm ../fft3d_wrap.cpp
  rm ../remap.cpp
  rm ../remap_wrap.cpp

  rm ../ewald.h
  rm ../pppm.h
  rm ../pppm_tip4p.h
  rm ../pair_born_coul_long.h
  rm ../pair_buck_coul_long.h
  rm ../pair_coul_long.h
  rm ../pair_lj_cut_coul_long.h
  rm ../pair_lj_cut_coul_long_tip4p.h
  rm ../pair_lj_charmm_coul_long.h
  rm ../fft3d.h
  rm ../fft3d_wrap.h
  rm ../remap.h
  rm ../remap_wrap.h

endif
