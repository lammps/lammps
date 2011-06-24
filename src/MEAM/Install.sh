# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude MEAM info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*meam[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/meam |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/meam |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lmeam |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(meam_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(meam_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(meam_SYSPATH) |' ../Makefile.package
  fi

  cp pair_meam.cpp ..

  cp pair_meam.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*meam[^ \t]* //' ../Makefile.package
  fi

  rm -f ../pair_meam.cpp

  rm -f ../pair_meam.h

fi
