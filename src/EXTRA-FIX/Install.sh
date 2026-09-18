# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0 || (test -n "$2" && test ! -e ../$2)) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    cp $1 ..
    if (test $mode = 2) then
      echo "  updating src/$1"
    fi
  fi
}

# Fix uvt requires the combined temperature compute from EXTRA-COMPUTE.

for file in *.cpp *.h; do
  case "$file" in
    fix_uvt.cpp|fix_uvt.h) action $file compute_temp_uvt.h ;;
    *) test -f ${file} && action $file ;;
  esac
done
