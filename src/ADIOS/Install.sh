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

  CONFIGSCRIPT=none
  if ( type adios2-config > /dev/null 2>&1 ) then
    CONFIGSCRIPT=adios2-config
  elif ( ! test -z "$ADIOS2_DIR" ) then
    if ( type $ADIOS2_DIR/bin/adios2-config > /dev/null 2>&1 ) then
      CONFIGSCRIPT=$ADIOS2_DIR/bin/adios2-config
    else
      echo "ERROR: ADIOS2_DIR environment variable is set but" \
           "\$ADIOS2_DIR/bin/adios2-config does not exist"
    fi
  elif ( ! test -z "$ADIOS_DIR" ) then
    if ( type $ADIOS_DIR/bin/adios2-config > /dev/null 2>&1 ) then
      CONFIGSCRIPT=$ADIOS_DIR/bin/adios2-config
    else
      echo "ERROR: ADIOS_DIR environment variable is set but" \
           "\$ADIOS_DIR/bin/adios2-config does not exist"
    fi
  else
    echo "ERROR: ADIOS2_DIR environment variable must point to ADIOS 2.x" \
         "installation directory or adios2-config should be in PATH"
  fi

  if [ "$CONFIGSCRIPT" != "none" ]; then
    ADIOS2_INC=`$CONFIGSCRIPT --cxx-flags`
    ADIOS2_LIB=`$CONFIGSCRIPT --cxx-libs`

    echo "adios_SYSINC=${ADIOS2_INC}
adios_SYSLIB=${ADIOS2_LIB}
" > Makefile.lammps


    if (test -e ../Makefile.package) then
      sed -i -e 's/[^ \t]*adios[^ \t]* //g' ../Makefile.package
      sed -i -e '/^adios_SYS.*$/d' ../Makefile.package
      sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(adios_SYSINC) |' ../Makefile.package
      sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(adios_SYSLIB) |' ../Makefile.package
    fi

    if (test -e ../Makefile.package.settings) then
      sed -i -e '/^[ \t]*include.*ADIOS.*$/d' ../Makefile.package.settings
      # multiline form needed for BSD sed on Macs
      sed -i -e '4 i \
include ../ADIOS/Makefile.lammps
' ../Makefile.package.settings
    fi
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*adios[^ \t]* //g' ../Makefile.package
    sed -i -e '/^adios_SYS.*$/d' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*ADIOS.*$/d' ../Makefile.package.settings
  fi

  rm -f Makefile.lammps

fi
