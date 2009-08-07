# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude MEAM library

if ($1 == 1) then

  cp style_meam.h ..

  cp pair_meam.cpp ..

  cp pair_meam.h ..

  sed -i 's/\S*meam //' ../Makefile.package
  sed -i 's|^PKGINC =\s*|&-I../../lib/meam |' ../Makefile.package
  sed -i 's|^PKGPATH =\s*|&-L../../lib/meam |' ../Makefile.package
  sed -i 's|^PKGLIB =\s*|&-lmeam |' ../Makefile.package

else if ($1 == 0) then

  rm ../style_meam.h
  touch ../style_meam.h

  rm ../pair_meam.cpp

  rm ../pair_meam.h

  sed -i 's/\S*meam //' ../Makefile.package

endif
