# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kim[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(kim_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(kim_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(kim_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*kim.*$/d' ../Makefile.package.settings
    sed -i '4 i include ..\/..\/lib\/kim\/Makefile.lammps' ../Makefile.package.settings
  fi

  cp pair_kim.cpp ..
  cp pair_kim.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*kim[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*kim.*$/d' ../Makefile.package.settings
  fi

  rm -f ../pair_kim.cpp
  rm -f ../pair_kim.h

fi
