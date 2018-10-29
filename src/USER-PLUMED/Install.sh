# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude ATC info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's|^PKG_LIB =[ \t]*|& -lplumedWrapper -ldl |' ../Makefile.package
    if ( ! test -e ../../lib/plumed/liblink/plumed/src/lib/Plumed.inc.static ) then
        sed -i -e 's|^PKG_SYSINC =[ \t]*|& -D__PLUMED_HAS_DLOPEN |' ../Makefile.package
    fi
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/plumed/includelink |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/plumed/liblink |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    # This is for statically linking plumed2
    if ( test -e ../../lib/plumed/liblink/src/lib/Plumed.inc.static ) then  
       fname=../../lib/plumed/liblink/src/lib/Plumed.inc.static 
       sed -i -e '4 i \
include '$fname'
' ../Makefile.package.settings
    # This is for linking plumed2 as a runtime library  -- this is the default behavior
    else 
    # multiline form needed for BSD sed on Macs
       sed -i -e '4 i \
PLUMED_LOAD=-ldl
' ../Makefile.package.settings
    fi
  fi

  cp fix_plumed.cpp ..
  cp fix_plumed.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*-lplumedWrapper -ldl[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*-D__PLUMED_HAS_DLOPEN[^ \t]* //' ../Makefile.package
    sed -i -e 's|[^ \t]*-I../../lib/plumed/includelink[^ \t]* ||' ../Makefile.package
    sed -i -e 's|[^ \t]*-L../../lib/plumed/liblink[^ \t]* ||' ../Makefile.package
  fi

  rm -f ../fix_plumed.cpp
  rm -f ../fix_plumed.h
fi
