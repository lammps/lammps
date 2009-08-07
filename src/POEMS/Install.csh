# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude POEMS library

if ($1 == 1) then

  cp style_poems.h ..

  cp fix_poems.cpp ..

  cp fix_poems.h ..

  sed -i 's/\S*poems //' ../Makefile.package
  sed -i 's|^PKGINC =\s*|&-I../../lib/poems |' ../Makefile.package
  sed -i 's|^PKGPATH =\s*|&-L../../lib/poems |' ../Makefile.package
  sed -i 's|^PKGLIB =\s*|&-lpoems |' ../Makefile.package

else if ($1 == 0) then

  rm ../style_poems.h
  touch ../style_poems.h

  rm ../fix_poems.cpp

  rm ../fix_poems.h

  sed -i 's/\S*poems //' ../Makefile.package

endif
