# Install/unInstall package files in LAMMPS

# step 1: process all *_omp.cpp and *_omp.h files.
# do not install child files if parent does not exist

for file in *_omp.cpp *_omp.h ; do
  ofile=`echo $file | sed  -e 's,\(.*\)_omp\.h,\1.h,' -e 's,\(.*\)_omp\.cpp,\1.cpp,'`
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

# step 2: handle cases and tasks not handled in step 1.
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_USER_OMP |' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_OMP switch

  touch ../accelerator_omp.h

  cp thr_data.h ..
  cp thr_data.cpp ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_OMP switch

  touch ../accelerator_omp.h

  rm -f ../thr_data.h
  rm -f ../thr_data.cpp

fi
