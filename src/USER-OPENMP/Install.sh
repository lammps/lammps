# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

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

  rm ../pair_omp.cpp
  rm ../pair_lj_cut_omp.cpp
  rm ../pair_gauss_cut_omp.cpp

  rm ../pair_omp.h
  rm ../pair_lj_cut_omp.h
  rm ../pair_gauss_cut_omp.h

  rm -f ../pair_lj_charmm_coul_long_omp.cpp
  rm -f ../pair_lj_charmm_coul_long_omp.h

fi
