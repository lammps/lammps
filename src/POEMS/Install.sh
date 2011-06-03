# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude POEMS info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*poems //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/poems |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/poems |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lpoems |' ../Makefile.package
  fi

  cp fix_poems.cpp ..

  cp fix_poems.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*poems //' ../Makefile.package
  fi

  rm -f ../fix_poems.cpp

  rm -f ../fix_poems.h

fi
