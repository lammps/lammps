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

# list of files with optional dependencies

action compute_temp_cs.cpp
action compute_temp_cs.h

action pair_born_coul_long_cs.cpp   pair_born_coul_long.cpp
action pair_born_coul_dsf_cs.cpp    pair_born_coul_dsf.cpp
action pair_buck_coul_long_cs.cpp   pair_buck_coul_long.cpp
action pair_born_coul_long_cs.h     pair_born_coul_long.h
action pair_born_coul_dsf_cs.h      pair_born_coul_dsf.h
action pair_buck_coul_long_cs.h     pair_buck_coul_long.h
action pair_coul_long_cs.cpp        pair_coul_long.cpp
action pair_coul_long_cs.h          pair_coul_long.h
action pair_lj_cut_coul_long_cs.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_cs.h   pair_lj_cut_coul_long.h
action pair_lj_class2_coul_long_cs.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_cs.h pair_lj_class2_coul_long.h

action pair_born_coul_wolf_cs.cpp   pair_born_coul_wolf.cpp
action pair_born_coul_wolf_cs.h     pair_born_coul_wolf.h
action pair_coul_wolf_cs.cpp        pair_coul_wolf.cpp
action pair_coul_wolf_cs.h          pair_coul_wolf.h
