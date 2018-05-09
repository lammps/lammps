# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# this is default Install.sh for all packages
# if package has an auxiliary library or a file with a dependency,
# then package dir has its own customized Install.sh

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

# all package files with dependencies

action fix_qeq.cpp
action fix_qeq.h
action fix_qeq_dynamic.cpp
action fix_qeq_dynamic.h
action fix_qeq_fire.cpp pair_comb.h
action fix_qeq_fire.h pair_comb.h
action fix_qeq_point.cpp
action fix_qeq_point.h
action fix_qeq_shielded.cpp
action fix_qeq_shielded.h
action fix_qeq_slater.cpp
action fix_qeq_slater.h
