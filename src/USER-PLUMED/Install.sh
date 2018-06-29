# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude ATC info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's|^PKG_LIB =[ \t]*|& $(PLUMED_LOAD) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/Plumed.inc
' ../Makefile.package.settings
  fi

  cp fix_plumed.cpp ..
  cp fix_plumed.h ..
  cp Plumed.h ..
  cp Plumed.cpp ..
  cp Plumed.inc ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]* \$(PLUMED_LOAD)[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*Plumed\.inc.*$/d' ../Makefile.package.settings
  fi

  rm -f ../fix_plumed.cpp
  rm -f ../fix_plumed.h
  rm -f ../Plumed.h
  rm -f ../Plumed.cpp
  rm -f ../Plumed.inc

fi
