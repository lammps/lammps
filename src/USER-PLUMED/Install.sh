# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude ATC info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's|^PKG_LIB =[ \t]*|& $(PLUMED_LOAD) |' ../Makefile.package
    if ( ! test -e ../../lib/plumed/plumed2*/src/lib/Plumed.inc.static ) then
        sed -i -e 's|^PKG_SYSINC =[ \t]*|& -D__PLUMED_HAS_DLOPEN |' ../Makefile.package
    fi
  fi

  if (test -e ../Makefile.package.settings) then
    # This is for statically linking plumed2
    if ( test -e ../../lib/plumed/plumed2*/src/lib/Plumed.inc.static ) then  
       fname=`ls ../../lib/plumed/plumed2*/src/lib/Plumed.inc.static` 
       sed -i -e '4 i \
include '$fname'
' ../Makefile.package.settings
       dname=`ls ../../lib/plumed/plumed2*/src/wrapper/Plumed.h`
       ln -s USER-PLUMED/$dname ../Plumed.h
    # This is for linking plumed2 as a runtime library  -- this is the default behavior
    else 
    # multiline form needed for BSD sed on Macs
       sed -i -e '4 i \
PLUMED_LOAD=-ldl
' ../Makefile.package.settings
       cp Plumed.h ..
       cp Plumed.cpp ..
    fi
  fi

  cp fix_plumed.cpp ..
  cp fix_plumed.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]* \$(PLUMED_LOAD)[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]* -D__PLUMED_HAS_DLOPEN[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    fname=`ls ../../lib/plumed/plumed2*/src/lib/Plumed.inc.static`
    sed -i -e "\:include $fname: d" ../Makefile.package.settings 
    sed -i -e '/PLUMED_LOAD=-ldl/d' ../Makefile.package.settings
  fi

  rm -f ../fix_plumed.cpp
  rm -f ../fix_plumed.h
  rm -f ../Plumed.h
  if ( test -e ../Plumed.h ) then
       rm -f ../Plumed.cpp
  fi

fi
