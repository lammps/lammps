#/bin/sh
# Update package files in LAMMPS
# copy package file to src if it doesn't exists or is different
# do not copy gayberne files if non-GPU version does not exist
for file in *_omp.cpp *_omp.h; do
  # let us see if the "rain man" can count the toothpicks...
  ofile=`echo $file | sed -e s,\\\\\\(.\\*\\\\\\)_omp\\\\.\\\\\\(h\\\\\\|cpp\\\\\\),\\\\1.\\\\2,`
  if (test $file = "thr_omp.h") || (test $file = "thr_omp.cpp") then
    :  # do check for those files.
  elif (test ! -e ../$ofile) then
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

