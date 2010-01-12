# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude POEMS library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*poems //' ../Makefile.package
  sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/poems |' ../Makefile.package
  sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/poems |' ../Makefile.package
  sed -i -e 's|^PKG_LIB =[ \t]*|&-lpoems |' ../Makefile.package

  cp fix_poems.cpp ..

  cp fix_poems.h ..

else if ($1 == 0) then

  sed -i -e 's/[^ \t]*poems //' ../Makefile.package

  rm ../fix_poems.cpp

  rm ../fix_poems.h

endif
