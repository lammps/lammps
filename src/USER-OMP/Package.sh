# Update package files in LAMMPS
# copy package file to src if it doesn't exists or is different
# do not copy OpenMP style files, if a non-OpenMP version does 
# not exist. Do remove OpenMP style files that have no matching
# non-OpenMP version installed, e.g. after a package has been
# removed
for file in *_omp.cpp *_omp.h pppm_proxy.h pppm_proxy.cpp thr_data.h thr_data.cpp; do
  # let us see if the "rain man" can count the toothpicks...
  ofile=`echo $file | sed -e s,\\\\\\(.\\*\\\\\\)_omp\\\\.\\\\\\(h\\\\\\|cpp\\\\\\),\\\\1.\\\\2,`
  if (test $file = "thr_omp.h") || (test $file = "thr_omp.cpp") \
      || (test $file = "thr_data.h") || (test $file = "thr_data.cpp") then
    if (test ! -e ../$file) then
      echo "  creating src/$file"
      cp $file ..
    elif ! cmp -s $file ../$file ; then
      echo "  updating src/$file"
      cp $file ..
    fi
  elif (test $file = "pppm_proxy.h") || (test $file = "pppm_proxy.cpp") \
      || (test $file = "pair_lj_charmm_coul_pppm_omp.h") \
      || (test $file = "pair_lj_charmm_coul_pppm_omp.cpp") then
    if (test ! -e ../pair_lj_charmm_coul_long.h ) then
      echo "  removing src/$file"
      rm -f ../$file
    else
      if (test ! -e ../$file) then
        echo "  creating src/$file"
        cp $file ..
      elif ! cmp -s $file ../$file ; then
        echo "  updating src/$file"
        cp $file ..
      fi
    fi
    continue
  elif (test ! -e ../$ofile) then
    if (test -e ../$file) then
      echo "  removing src/$file"
      rm -f ../$file
    fi
    continue
  else
    if (test ! -e ../$file) then
      echo "  creating src/$file"
      cp $file ..
    elif ! cmp -s $file ../$file ; then
      echo "  updating src/$file"
      cp $file ..
    fi
  fi
done

for file in thr_data.h thr_data.cpp; do
  if (test ! -e ../$file) then
    echo "  creating src/$file"
    cp $file ..
  elif ! cmp -s $file ../$file ; then
    echo "  updating src/$file"
    cp $file ..
  fi
done

