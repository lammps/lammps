# Install/unInstall package files in LAMMPS

if (test $1 = 1) then

  echo "================================================================================"
  echo "USER-SELM Package: Fluctuating Hydrodynamics Package"
  echo "--------------------------------------------------------------------------------"

  ## PJA: can be helpful in some cases to use symbolic links.
  #echo "Symbolic link copying files for USER-SELM package."    
  #basePath=$(pwd -P)
  #echo "Base Path = $basePath"
  #ln -sf $basePath/*.h ../
  #ln -sf $basePath/*.cpp ../

  echo "Copying files for the package into the source directory."
  #cp -p $PWD/*.h ../
  #cp -p $PWD/*.cpp ../
  find $PWD -name '*.h' -exec cp "{}" ../ >& selm_cp_h.log \;
  find $PWD -name '*.cpp' -exec cp "{}" ../ >& selm_cp_cpp.log \;

  echo " "
  echo "NOTE: This USER-SELM package version Should be run in serial mode."
  echo "Uses the head node for handling the fluid mechanics.  "
  echo " "
  echo "For more information and examples see "
  echo "http://mango-selm.org"
  echo " "
  echo "================================================================================"

elif (test $1 = 0) then

  rm ../fix_selm.cpp
  rm ../fix_selm.h

  echo "  "
#  echo "WARNING: (PJA) List of files to remove is not yet implemented."

fi


# test is MOLECULE package installed already
if (test $1 = 1) then
  if (test ! -e ../angle_harmonic.cpp) then
    echo "Must install MOLECULE package with USER-SELM, for example by command 'make yes-molecule'."
    exit 1
  fi

fi

# setup library settings
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*selm[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/selm |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/selm |' ../Makefile.package

#   sed -i -e 's|^PKG_INC =[ \t]*|&-I../../lib/selm |' ../Makefile.package
#   sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/selm$(LIBSOBJDIR) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lselm |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(user-selm_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-selm_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-selm_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*selm.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/selm\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*selm[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*selm.*$/d' ../Makefile.package.settings
  fi

fi

