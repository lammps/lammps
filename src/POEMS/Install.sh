# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude POEMS info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*poems[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/poems |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/poems |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lpoems |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*poems.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/poems\/Makefile.lammps\
' ../Makefile.package.settings
  fi

  cp fix_poems.cpp ..

  cp fix_poems.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*poems[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*poems.*$/d' ../Makefile.package.settings
  fi

  rm -f ../fix_poems.cpp

  rm -f ../fix_poems.h

fi
