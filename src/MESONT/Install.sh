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

# list of files with optional dependencies
action angle_mesocnt.cpp
action angle_mesocnt.h
action bond_mesocnt.cpp bond_harmonic.cpp
action bond_mesocnt.h bond_harmonic.h
action compute_mesont.cpp
action compute_mesont.h
action pair_mesocnt.cpp
action pair_mesocnt.h
action pair_mesocnt_viscous.cpp
action pair_mesocnt_viscous.h

action export_mesont.h
action atom_vec_mesont.cpp
action atom_vec_mesont.h
action pair_mesont_tpm.cpp
action pair_mesont_tpm.h

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*mesont[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/mesont |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/mesont |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lmesont |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(mesont_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(mesont_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(mesont_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*mesont.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/mesont\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*mesont[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*mesont.*$/d' ../Makefile.package.settings
  fi

fi
