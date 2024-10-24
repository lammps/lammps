# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

cat <<EOF
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING

  Support for building the COLVARS package with the legacy build system using GNU
  make will be removed in Summer 2025.  Please switch to using CMake to build
  LAMMPS as soon as possible and report any problems to developers@lammps.org
  or post a bug report issue at https://github.com/lammps/lammps/issues

WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING-WARNING
EOF

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
    sed -i -e 's/[^ \t]*colvars[^ \t]* //g' ../Makefile.package
    if (test ! -e ../pair_lepton.cpp) then
      sed -i -e 's/[^ \t]*lepton[^ \t]* //g' ../Makefile.package
    fi
    sed -i -e 's|^PKG_INC =[ \t]*|&-I..\/..\/lib\/colvars -I..\/..\/lib\/lepton\/include -I..\/..\/lib\/lepton |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/colvars$(LIBOBJDIR) -L..\/..\/lib\/lepton$(LIBOBJDIR) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lcolvars -llepton |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(colvars_SYSINC) $(lepton_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(colvars_SYSLIB) $(lepton_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(colvars_SYSPATH) $(lepton_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*colvars.*$/d' ../Makefile.package.settings
    if (test ! -e ../pair_lepton.cpp) then
      sed -i -e '/^[ \t]*include.*lepton.*$/d' ../Makefile.package.settings
    fi
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/colvars\/Makefile.lammps
' ../Makefile.package.settings

    sed -i -e '/^[ \t]*include.*lepton.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/lepton\/Makefile.lammps
' ../Makefile.package.settings

  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*colvars[^ \t]* //g' ../Makefile.package
    if (test ! -e ../pair_lepton.cpp) then
      sed -i -e 's/[^ \t]*lepton[^ \t]* //g' ../Makefile.package
    fi
  fi
  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*colvars.*$/d' ../Makefile.package.settings
    if (test ! -e ../pair_lepton.cpp) then
      sed -i -e '/^[ \t]*include.*lepton.*$/d' ../Makefile.package.settings
    fi
  fi
fi
