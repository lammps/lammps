# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude CUDA info
# do not install child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*colvars[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I..\/..\/lib\/colvars |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/colvars |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lcolvars |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(colvars_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(colvars_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(colvars_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*colvars.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/colvars\/Makefile.lammps\
' ../Makefile.package.settings

  fi

  cp colvarproxy_lammps.h ..
  cp colvarproxy_lammps.cpp ..
  cp fix_colvars.h ..
  cp fix_colvars.cpp ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*colvars[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*colvars.*$/d' ../Makefile.package.settings
  fi

  rm -f ../colvarproxy_lammps.h
  rm -f ../colvarproxy_lammps.cpp
  rm -f ../fix_colvars.h
  rm -f ../fix_colvars.cpp
fi
