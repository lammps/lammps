# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  cp -p pair_gauss_cut.cpp ..
  cp -p pair_coul_diel.cpp ..

  cp -p pair_gauss_cut.h ..
  cp -p pair_coul_diel.h ..

#  if (test -e ../pair_lj_charmm_coul_long.cpp) then
#      cp -p pair_lj_charmm_coul_long_omp.cpp ..
#      cp -p pair_lj_charmm_coul_long_omp.h ..
#  fi

elif (test $1 = 0) then

  rm ../pair_gauss_cut.cpp
  rm ../pair_coul_diel.cpp

  rm ../pair_gauss_cut.h
  rm ../pair_coul_diel.h

#  rm -f ../pair_lj_charmm_coul_long_omp.cpp
#  rm -f ../pair_lj_charmm_coul_long_omp.h

fi
