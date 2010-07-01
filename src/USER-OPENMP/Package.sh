# Update package files in LAMMPS
# cp package file to src if doesn't exist or is different
# do not copy any package files if non-OMP version does not exist

for file in *_omp.cpp *_omp.h; do
   # let us see if the "rain man" can count the toothpicks... 
   ofile=`echo $file | sed -e s,\\\\\\(.\\*\\\\\\)_omp\\\\.\\\\\\(h\\\\\\|cpp\\\\\\),\\\\1.\\\\2,`
   if (test ! -e ../$ofile) then
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
