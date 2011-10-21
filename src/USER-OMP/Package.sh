# Update package files in LAMMPS
# cp package file to src if doesn't exist or is different
# do not copy certain files if non-OMP versions do not exist
# do remove OpenMP style files that have no matching
#   non-OpenMP version installed, e.g. after a package has been removed

for file in *_omp.cpp *_omp.h; do
  # let us see if the "rain man" can count the toothpicks...
  ofile=`echo $file | sed -e s,\\\\\\(.\\*\\\\\\)_omp\\\\.\\\\\\(h\\\\\\|cpp\\\\\\),\\\\1.\\\\2,`
  if (test $file = "thr_omp.h") || (test $file = "thr_omp.cpp") then
    :  # always check for those files.
  elif (test ! -e ../$ofile) then
    if (test -e ../$file) then
      echo "  removing src/$file"
      rm -f ../$file
    fi
    continue
  fi

  if (test ! -e ../$file) then
    echo "  creating src/$file"
    cp $file ..
  elif ! cmp -s $file ../$file ; then
    echo "  updating src/$file"
    cp $file ..
  fi
done

