# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

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

for file in *.cpp *.h; do
  action $file
done

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

#  if (test -z "$ADIOS_DIR") then
#    if command -v adios_config; then 
#        ADIOS_DIR=`adios_config -d`
#    else
#      echo "ERROR: ADIOS_DIR environment variable needs to point to ADIOS" \
#           " installation directory or adios_config should be in PATH"
#    fi
#  fi
#  ADIOS_INC=-I${ADIOS_DIR}/include
#  ADIOS_LIB=`${ADIOS_DIR}/bin/adios_config -l`
#  
#  echo "adios_SYSINC=${ADIOS_INC}
#adios_SYSLIB=${ADIOS_LIB}
#adios_SYSPATH=${ADIOS_DIR}
#" > ../Makefile.adios 


  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*adios[^ \t]* //g' ../Makefile.package
    sed -i -e 's/-DLMP_ADIOS //g' ../Makefile.package
    sed -i -e '/^adios_SYS.*$/d' ../Makefile.package
#    sed -i -e '4 i \
#adios_SYSINC='"${ADIOS_INC}"'
#' ../Makefile.package
#    sed -i -e '5 i \
#adios_SYSLIB='"${ADIOS_LIB}"'
#' ../Makefile.package
#    sed -i -e '6 i \
#adios_SYSPATH='"${ADIOS_DIR}"'
#' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_ADIOS |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(adios_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(adios_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(adios_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*adios.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/adios\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*adios[^ \t]* //g' ../Makefile.package
    sed -i -e 's/-DLMP_ADIOS //g' ../Makefile.package
    sed -i -e '/^adios_SYS.*$/d' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*adios.*$/d' ../Makefile.package.settings
  fi

fi
