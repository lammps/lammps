# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude MEAM library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*meam //' ../Makefile.package
  sed -i -e 's/[^ \t]*meam_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/meam |' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/meam |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-lmeam |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(meam_SYSPATH) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(meam_SYSLIB) |' ../Makefile.package

  cp style_meam.h ..

  cp pair_meam.cpp ..

  cp pair_meam.h ..

else if ($1 == 0) then

  sed -i -e 's/[^ \t]*meam //' ../Makefile.package
  sed -i -e 's/[^ \t]*meam_[^ \t]*) //' ../Makefile.package

  rm ../style_meam.h
  touch ../style_meam.h

  rm ../pair_meam.cpp

  rm ../pair_meam.h

endif
