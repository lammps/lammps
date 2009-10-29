# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude ATC library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*atc //' ../Makefile.package
  sed -i -e 's/[^ \t]*atc_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/atc |' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/atc |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-latc |' ../Makefile.package
  sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-atc_SYSPATH) |' ../Makefile.package
  sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-atc_SYSLIB) |' ../Makefile.package

  cp style_user_atc.h ..

  cp fix_atc.h ..
  cp fix_atc.cpp ..

else if ($1 == 0) then

  sed -i -e 's/[^ \t]*atc //' ../Makefile.package
  sed -i -e 's/[^ \t]*atc_[^ \t]*) //' ../Makefile.package

  rm ../style_user_atc.h
  touch ../style_user_atc.h

  rm ../fix_atc.h
  rm ../fix_atc.cpp

endif

