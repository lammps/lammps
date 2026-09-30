# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

action () {
  if (test $mode = 0) then
    rm -f ../$1
  fi
}

# force rebuild of files using python header

touch ../lmppython.h

# all package files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1 || test $1 = 2) then
  echo "The PYTHON package no longer supports the legacy build system. Please build LAMMPS with CMake instead."
  exit 1

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*python[^ \t]* //' ../Makefile.package
    sed -i -e 's/[^ \t]*PYTHON[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*python.*$/d' ../Makefile.package.settings
  fi

fi
