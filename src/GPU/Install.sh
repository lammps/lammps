# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude GPU library
# do not copy gayberne files if non-GPU version does not exist

if (test $1 = 1) then

  sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
  sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package

  if (test -e ../pair_gayberne.cpp) then
    cp pair_gayberne_gpu.cpp ..
    cp pair_gayberne_gpu.h ..
  fi

  cp pair_lj_cut_gpu.cpp ..
  cp pair_lj_cut_gpu.h ..

elif (test $1 = 0) then

  sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
  sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package

  rm ../pair_gayberne_gpu.cpp
  rm ../pair_lj_cut_gpu.cpp

  rm ../pair_gayberne_gpu.h
  rm ../pair_lj_cut_gpu.h

fi
