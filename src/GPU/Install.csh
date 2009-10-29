# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude GPU library
# do not copy gayberne files if non-GPU version does not exist

if ($1 == 1) then

  sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
  sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package

  cp style_gpu.h tmp.h
  if (! -e ../pair_gayberne.cpp) then
    grep -v gayberne tmp.h > tmp1.h
    mv tmp1.h tmp.h
  endif
  mv tmp.h ../style_gpu.h

  if (-e ../pair_gayberne.cpp) then
    cp pair_gayberne_gpu.cpp ..
    cp pair_gayberne_gpu.h ..
  endif

  cp pair_lj_cut_gpu.cpp ..
  cp pair_lj_cut_gpu.h ..

else if ($1 == 0) then

  sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
  sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package

  rm ../style_gpu.h
  touch ../style_gpu.h

  rm ../pair_gayberne_gpu.cpp
  rm ../pair_lj_cut_gpu.cpp

  rm ../pair_gayberne_gpu.h
  rm ../pair_lj_cut_gpu.h

endif
