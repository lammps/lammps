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

# only uninstall
for file in *.cpp *.h; do
    action ${file}
done

if (test $1 = 1 || test $1 = 2) then
  echo "The APIP package no longer supports the legacy build system. Please build LAMMPS with CMake instead."
  exit 1
fi
