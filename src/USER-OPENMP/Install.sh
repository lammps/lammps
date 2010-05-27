# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  #Added 05/25/10
  cp pair_lj96_cut_omp.cpp ..
  cp pair_lj96_cut_omp.h ..
  cp pair_lj_smooth_omp.cpp ..
  cp pair_lj_smooth_omp.h ..
  cp pair_buck_omp.cpp ..
  cp pair_buck_omp.h ..

  #Added 05/26/10
  cp pair_morse_omp.cpp ..
  cp pair_morse_omp.h ..
  cp pair_born_coul_long_omp.cpp ..
  cp pair_born_coul_long_omp.h ..
  cp pair_soft_omp.cpp ..
  cp pair_soft_omp.h ..
  cp pair_yukawa_omp.cpp ..
  cp pair_yukawa_omp.h ..

  #Added 05/27/10
  cp pair_lj_cut_coul_cut_omp.cpp ..
  cp pair_lj_cut_coul_cut_omp.h ..
  cp pair_lj_cut_coul_debye_omp.cpp ..
  cp pair_lj_cut_coul_debye_omp.h ..
  cp pair_lj_charmm_coul_charmm_omp.cpp ..
  cp pair_lj_charmm_coul_charmm_omp.h ..
  cp pair_lj_charmm_coul_charmm_implicit_omp.cpp ..
  cp pair_lj_charmm_coul_charmm_implicit_omp.h ..

  cp pair_omp.cpp ..
  cp pair_lj_cut_omp.cpp ..
  cp pair_gauss_cut_omp.cpp ..

  cp pair_omp.h ..
  cp pair_lj_cut_omp.h ..
  cp pair_gauss_cut_omp.h ..

  if (test -e ../pair_lj_charmm_coul_long.cpp) then
      cp -p pair_lj_charmm_coul_long_omp.cpp ..
      cp -p pair_lj_charmm_coul_long_omp.h ..
  fi

elif (test $1 = 0) then

  #Added 05/25/10
  rm ../pair_lj96_cut_omp.cpp
  rm ../pair_lj96_cut_omp.h
  rm ../pair_lj_smooth_omp.cpp
  rm ../pair_lj_smooth_omp.h
  rm ../pair_buck_omp.cpp 
  rm ../pair_buck_omp.h

  #Added 05/26/10
  rm ../pair_morse_omp.cpp 
  rm ../pair_morse_omp.h
  rm ../pair_born_coul_long_omp.cpp
  rm ../pair_born_coul_long_omp.h
  rm ../pair_soft_omp.cpp 
  rm ../pair_soft_omp.h
  rm ../pair_yukawa_omp.cpp
  rm ../pair_yukawa_omp.h

  #Added 05/27/10
  rm ../pair_lj_cut_coul_cut_omp.cpp 
  rm ../pair_lj_cut_coul_cut_omp.h
  rm ../pair_lj_cut_coul_debye_omp.cpp
  rm ../pair_lj_cut_coul_debye_omp.h
  rm ../pair_omp.cpp
  rm ../pair_lj_cut_omp.cpp
  rm ../pair_gauss_cut_omp.cpp
  rm ../pair_lj_charmm_coul_charmm_omp.cpp 
  rm ../pair_lj_charmm_coul_charmm_omp.h 
  rm ../pair_lj_charmm_coul_charmm_implicit_omp.cpp 
  rm ../pair_lj_charmm_coul_charmm_implicit_omp.h 

  rm ../pair_omp.h
  rm ../pair_lj_cut_omp.h
  rm ../pair_gauss_cut_omp.h

  rm -f ../pair_lj_charmm_coul_long_omp.cpp
  rm -f ../pair_lj_charmm_coul_long_omp.h

fi
