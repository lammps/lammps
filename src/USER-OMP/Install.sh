# Install/unInstall package files in LAMMPS
# do not install child files if parent does not exist

if (test $1 = 1) then

#  if (test -e ../pair_lj_cut_coul_long.cpp) then
#    cp pair_lj_cut_coul_long_omp.cpp ..
#    cp pair_lj_cut_coul_long_omp.h ..
#  fi

  cp pair_lj_cut_omp.cpp ..

  cp thr_omp.cpp ..

  cp pair_lj_cut_omp.h ..

  cp thr_omp.h ..

elif (test $1 = 0) then

#  rm -f ../pair_lj_cut_coul_long_omp.cpp
  rm -f ../pair_lj_cut_omp.cpp

  rm -f ../thr_omp.cpp

#  rm -f ../pair_lj_cut_coul_long_omp.h
  rm -f ../pair_lj_cut_omp.h

  rm -f ../thr_omp.h

fi

