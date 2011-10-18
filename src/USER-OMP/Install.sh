# Install/unInstall package files in LAMMPS
# do not install child files if parent does not exist

for file in *_omp.cpp *_omp.h; do
    # let us see if the "rain man" can count the toothpicks...
   ofile=`echo $file | sed -e s,\\\\\\(.\\*\\\\\\)_omp\\\\.\\\\\\(h\\\\\\|cpp\\\\\\),\\\\1.\\\\2,`

  if (test $1 = 1) then

    if (test $file = "thr_omp.h") || (test $file = "thr_omp.cpp") then
      :  # always install those files.
    elif (test ! -e ../$ofile) then
      continue
    fi

    cp $file ..

  elif (test $1 = 0) then

    rm -f ../$file
  fi
done


if (test $1 = 1) then

  cp thr_data.h ..
  cp thr_data.cpp ..

  cp pppm_proxy.h ..
  cp pppm_proxy.cpp ..

  cp pair_lj_charmm_coul_pppm_omp.h ..
  cp pair_lj_charmm_coul_pppm_omp.cpp ..

elif (test $1 = 0) then

  rm -f ../thr_data.h
  rm -f ../thr_data.cpp

  rm -f ../pppm_proxy.h
  rm -f ../pppm_proxy.cpp

  rm -f ../pair_lj_charmm_coul_pppm_omp.h
  rm -f ../pair_lj_charmm_coul_pppm_omp.cpp
fi
