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

# USER-PHONON uses the parallel FFT wrapper used in PPPM,
# so we must require the KSPACE package to be installed.

if (test $1 = 1) then
  if (test ! -e ../fft3d_wrap.h) then
    echo "Must install KSPACE package with USER-PHONON"
    exit 1
  fi
fi

# list of files with optional dependcies

action fix_phonon.cpp fft3d_wrap.h
action fix_phonon.h fft3d_wrap.h
action dynamical_matrix.cpp
action dynamical_matrix.h
