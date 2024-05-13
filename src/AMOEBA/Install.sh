# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

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

# pair style amoeba calls KSPACE functions and requires FFT grid.

if (test $1 = 1) then
  if (test ! -e ../pppm.cpp) then
    echo "Must install KSPACE package with AMOEBA package"
    exit 1
  fi
fi

for file in *.cpp *.h; do
    action ${file}
done
