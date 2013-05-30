# Install/unInstall/Update package files in LAMMPS
mode=$1

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (test -z "$2" || test -e ../$2) then
    if (test $mode = 1) then
      cp $1 ..
    elif (test $mode = 2) then
      if (! cmp -s $1 ../$1) then
        echo "updating src/$1"
        cp $1 ..
      fi
    fi
  fi
}

# step 1: process all *_omp.cpp and *_omp.h files.
# do not install child files if parent does not exist

for file in *_omp.cpp; do
  test $file = thr_omp.cpp && continue
  dep=${file%_omp.cpp}.cpp
  action $file $dep
done

for file in *_omp.h; do
  test $file = thr_omp.h && continue
  dep=${file%_omp.h}.h
  action $file $dep
done

action thr_omp.h
action thr_omp.cpp
action thr_data.h
action thr_data.cpp


# step 2: handle cases and tasks not handled in step 1.
if (test $mode = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_USER_OMP |' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_OMP switch

  touch ../accelerator_omp.h

elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_OMP switch

  touch ../accelerator_omp.h

fi
