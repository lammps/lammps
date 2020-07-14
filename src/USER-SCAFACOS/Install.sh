# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/`.*scafacos.*` //' ../Makefile.package
    sed -i -e 's/[^ \t]*scafacos[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/scafacos/includelink |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/scafacos/liblink |' ../Makefile.package
    sed -i -e 's%^PKG_LIB =[ \t]*%&`grep Libs: ../../lib/scafacos/liblink/pkgconfig/scafacos.pc | cut -d " " -f 2-` %' ../Makefile.package
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/`.*scafacos.*` //' ../Makefile.package
    sed -i -e 's/[^ \t]*scafacos[^ \t]* //' ../Makefile.package
  fi

fi
