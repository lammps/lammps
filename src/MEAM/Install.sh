# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude package info

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

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*meam.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/meam\/Makefile.lammps\
' ../Makefile.package.settings
  fi

  cp pair_meam.cpp ..

  cp pair_meam.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*meam[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*meam.*$/d' ../Makefile.package.settings
  fi

  rm -f ../pair_meam.cpp

  rm -f ../pair_meam.h

fi
