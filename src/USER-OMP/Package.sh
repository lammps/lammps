# Update package files in LAMMPS
# Copy package file to src if it doesn't exists or is different.
# But only copy the file, if a non-OpenMP version exists and
# remove OpenMP versions that have no matching serial file
# installed, e.g. after a package has been removed.
for file in *_omp.cpp *_omp.h ; do
  # these are special cases and handled below
  if (test $file = "thr_omp.h") || (test $file = "thr_omp.cpp") then
    continue
  fi
  # derive name of non-OpenMP version
  ofile=`echo $file | sed  -e 's,\(.*\)_omp\.h,\1.h,' -e 's,\(.*\)_omp\.cpp,\1.cpp,'`
  if (test ! -e ../$ofile) then
    if (test -e ../$file) then
      echo "  removing src/$file"
      rm -f ../$file
    fi
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

# special case for files not covered by the automatic script above
for file in thr_data.h thr_data.cpp thr_omp.h thr_omp.cpp; do
  if (test ! -e ../$file) then
    echo "  creating src/$file"
    cp $file ..
  elif ! cmp -s $file ../$file ; then
    echo "  updating src/$file"
    cp $file ..
  fi
done

