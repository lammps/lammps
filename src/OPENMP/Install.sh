# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
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
  test $file = reaxff_omp.h && continue
  dep=${file%_omp.h}.h
  action $file $dep
done

action reaxff_omp.h reaxff_api.h
action thr_omp.h
action thr_omp.cpp
action thr_data.h
action thr_data.cpp

# step 2: handle cases and tasks not handled in step 1

if (test $mode = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_OPENMP |' ../Makefile.package
  fi

  # need to delete a bunch of dependency files because they
  # indirectly depend on user_cuda.h

  for f in finish.d modify_cuda.d
  do \
    rm -f ../Obj_*/$f
  done

  # force rebuild of files with LMP_OPENMP switch

  touch ../accelerator_omp.h

elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*OMP[^ \t]* //' ../Makefile.package
  fi

  # need to delete a bunch of dependency files because they
  # indirectly depend on user_cuda.h

  for f in finish.d modify_cuda.d
  do \
    rm -f ../Obj_*/$f
  done

  # force rebuild of files with LMP_OPENMP switch

  touch ../accelerator_omp.h

fi
