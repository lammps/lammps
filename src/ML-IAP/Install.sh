# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

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

# enforce package dependency
if (test $1 = 1 || test $1 = 2) then
  if (test ! -e ../sna.h) then
     echo "Must install SNAP package to use ML-IAP package"
     exit 1
  fi
fi

# all package C++ files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# Install cython pyx file only if also Python is available
action mliap_model_python_couple.pyx python_impl.cpp

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then
  if (test "$(type cythonize 2> /dev/null)" != "" && test -e ../python_impl.cpp) then
    if (test -e ../Makefile.package) then
      sed -i -e 's|^PKG_INC =[ \t]*|&-DMLIAP_PYTHON |' ../Makefile.package
    fi
    if (test -e ../Makefile.package.settings) then
      sed -i -e '/^include.*python.*mliap_python.*$/d' ../Makefile.package.settings
      # multiline form needed for BSD sed on Macs
      sed -i -e '4 i \
include ..\/..\/lib\/python\/Makefile.mliap_python
' ../Makefile.package.settings
    fi
    cythonize -3 ../mliap_model_python_couple.pyx
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*-DMLIAP_PYTHON[^ \t]* //g' ../Makefile.package
  fi
  rm -f ../mliap_model_python_couple.cpp ../mliap_model_python_couple.h
  sed -i -e '/^include.*python.*mliap_python.*$/d' ../Makefile.package.settings

elif (test $1 = 2) then
  if (test "$(type cythonize 2> /dev/null)" != "" && test -e ../python_impl.cpp) then
    cythonize -3 ../mliap_model_python_couple.pyx
  else
    rm -f ../mliap_model_python_couple.cpp ../mliap_model_python_couple.h
  fi
fi
