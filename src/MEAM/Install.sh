# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude MEAM library

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*meam //' ../Makefile.package
    sed -i -e 's/[^ \t]*meam_[^ \t]*) //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/meam |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/meam |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lmeam |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(meam_SYSPATH) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(meam_SYSLIB) |' ../Makefile.package
  fi

  cp pair_meam.cpp ..

  cp pair_meam.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*meam //' ../Makefile.package
    sed -i -e 's/[^ \t]*meam_[^ \t]*) //' ../Makefile.package
  fi

  rm ../pair_meam.cpp

  rm ../pair_meam.h

fi
