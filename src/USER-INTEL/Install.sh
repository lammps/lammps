# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

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

# step 1: process all *_intel.cpp and *_intel.h files.
# do not install child files if parent does not exist

for file in *_intel.cpp; do
  test $file = thr_intel.cpp && continue
  dep=`echo $file | sed 's/neigh_full_intel/neigh_full/g' | \
      sed 's/_offload_intel//g' | sed 's/_intel//g'`
  action $file $dep
done

for file in *_intel.h; do
  test $file = thr_intel.h && continue
  dep=`echo $file | sed 's/_offload_intel//g' | sed 's/_intel//g'`
  action $file $dep
done

action intel_preprocess.h
action intel_buffers.h
action intel_buffers.cpp
action math_extra_intel.h

# step 2: handle cases and tasks not handled in step 1.

if (test $mode = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*INTEL[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_USER_INTEL |' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_INTEL switch

  touch ../accelerator_intel.h

elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*INTEL[^ \t]* //' ../Makefile.package
  fi

  # force rebuild of files with LMP_USER_INTEL switch

  touch ../accelerator_intel.h

fi

# step 3: map omp styles that are not in the intel package to intel suffix

#if (test $mode = 0) then
#
#  rm -f ../*ompinto_intel*
#
#else
#
#  echo "  The 'intel' suffix will use the USER-OMP package for all"
#  echo "  angle, bond, dihedral, kspace, and improper styles:"
#  stylelist="pair fix angle bond dihedral improper"
#  for header in $stylelist; do
#    HEADER=`echo $header | sed 's/\(.*\)/\U\1/'`
#    outfile=../$header"_ompinto_intel.h"
#    echo "    Creating $header style map: $outfile"
#    echo -n "// -- Header to map USER-OMP " > $outfile  
#    echo "styles to the intel suffix" >> $outfile
#    echo >> $outfile
#    echo "#ifdef "$HEADER"_CLASS" >> $outfile
#    grep -h 'Style(' ../$header*_omp.h | grep -v 'charmm/coul/long' | \
#	grep -v 'lj/cut' | grep -v 'gayberne' | \
#	sed 's/\/omp/\/intel/g' >> $outfile
#    echo "#endif" >> $outfile
#  done
#
#  header="kspace"
#  HEADER="KSPACE"
#  outfile=../$header"_ompinto_intel.h"
#  echo "    Creating $header style map: $outfile"
#  echo -n "// -- Header to map USER-OMP " > $outfile  
#  echo "styles to the intel suffix" >> $outfile
#  echo >> $outfile
#  echo "#ifdef "$HEADER"_CLASS" >> $outfile
#  grep -h 'KSpaceStyle(' ../*_omp.h | sed 's/\/omp/\/intel/g' >> $outfile
#  echo "#endif" >> $outfile
#
#fi
