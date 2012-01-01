# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp ewald.cpp ..
  cp pppm.cpp ..
  cp pppm_cg.cpp ..
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
  cp kissfft.h ..
  cp pppm.h ..
  cp pppm_cg.h ..
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

elif (test $1 = 0) then

  rm -f ../ewald.cpp
  rm -f ../pppm.cpp
  rm -f ../pppm_cg.cpp
  rm -f ../pppm_tip4p.cpp
  rm -f ../pair_born_coul_long.cpp
  rm -f ../pair_buck_coul_long.cpp
  rm -f ../pair_coul_long.cpp
  rm -f ../pair_lj_cut_coul_long.cpp
  rm -f ../pair_lj_cut_coul_long_tip4p.cpp
  rm -f ../pair_lj_charmm_coul_long.cpp
  rm -f ../fft3d.cpp
  rm -f ../fft3d_wrap.cpp
  rm -f ../remap.cpp
  rm -f ../remap_wrap.cpp

  rm -f ../ewald.h
  rm -f ../kissfft.h
  rm -f ../pppm.h
  rm -f ../pppm_cg.h
  rm -f ../pppm_tip4p.h
  rm -f ../pair_born_coul_long.h
  rm -f ../pair_buck_coul_long.h
  rm -f ../pair_coul_long.h
  rm -f ../pair_lj_cut_coul_long.h
  rm -f ../pair_lj_cut_coul_long_tip4p.h
  rm -f ../pair_lj_charmm_coul_long.h
  rm -f ../fft3d.h
  rm -f ../fft3d_wrap.h
  rm -f ../remap.h
  rm -f ../remap_wrap.h

fi
