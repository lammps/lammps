# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

action () {
  if (test "$mode" = 0); then
    rm -f ../$1
  fi
}

for file in *.cpp *.h; do
  test -f "$file" && action "$file"
done

if (test "$mode" = 1 || test "$mode" = 2); then
  echo "The QMMM-XTB package supports only the CMake build system."
  exit 1
fi
