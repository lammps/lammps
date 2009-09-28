# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude MEAM library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*meam //' ../Makefile.package
  sed -i -e 's/[^ \t]*meam_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKGINC =[ \t]*|&-I../../lib/meam |' ../Makefile.package
  sed -i -e 's|^PKGPATH =[ \t]*|&-L../../lib/meam |' ../Makefile.package
  sed -i -e 's|^PKGLIB =[ \t]*|&-lmeam |' ../Makefile.package
  sed -i -e 's|^PKGPATHSYS =[ \t]*|&$(meam_SYSLIBPATH) |' ../Makefile.package
  sed -i -e 's|^PKGLIBSYS =[ \t]*|&$(meam_SYSLIB) |' ../Makefile.package

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
