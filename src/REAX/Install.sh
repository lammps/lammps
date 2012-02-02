# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude REAX info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*reax[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/reax |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/reax |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lreax |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(reax_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(reax_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(reax_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*reax.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/reax\/Makefile.lammps\
' ../Makefile.package.settings
  fi

  cp pair_reax.cpp ..
  cp pair_reax.h ..
  cp pair_reax_fortran.h ..

  cp fix_reax_bonds.h ..
  cp fix_reax_bonds.cpp ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*reax[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*reax.*$/d' ../Makefile.package.settings
  fi

  rm -f ../pair_reax.cpp
  rm -f ../pair_reax.h
  rm -f ../pair_reax_fortran.h

  rm -f ../fix_reax_bonds.h
  rm -f ../fix_reax_bonds.cpp

fi
