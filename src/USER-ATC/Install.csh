# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude ATC library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*atc //' ../Makefile.package
  sed -i -e 's/[^ \t]*atc_[^ \t]*) //' ../Makefile.package
  sed -i -e 's|^PKGINC =[ \t]*|&-I../../lib/atc |' ../Makefile.package
  sed -i -e 's|^PKGPATH =[ \t]*|&-L../../lib/atc |' ../Makefile.package
  sed -i -e 's|^PKGLIB =[ \t]*|&-latc |' ../Makefile.package
  sed -i -e 's|^PKGPATHSYS =[ \t]*|&$(user-atc_SYSLIBPATH) |' ../Makefile.package
  sed -i -e 's|^PKGLIBSYS =[ \t]*|&$(user-atc_SYSLIB) |' ../Makefile.package

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

