# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude ATC info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*atc //' ../Makefile.package
    sed -i -e 's/[^ \t]*atc_[^ \t]*) //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/atc |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/atc |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-latc |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(user-atc_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-atc_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-atc_SYSPATH) |' ../Makefile.package
  fi

  cp fix_atc.h ..
  cp fix_atc.cpp ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*atc //' ../Makefile.package
    sed -i -e 's/[^ \t]*atc_[^ \t]*) //' ../Makefile.package
  fi

  rm -f ../fix_atc.h
  rm -f ../fix_atc.cpp

fi

