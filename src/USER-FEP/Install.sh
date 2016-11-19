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

# all package files with dependencies

action compute_fep.cpp
action compute_fep.h
action fix_adapt_fep.cpp
action fix_adapt_fep.h
action pair_coul_cut_soft.cpp
action pair_coul_cut_soft.h
action pair_coul_long_soft.cpp            pppm.cpp
action pair_coul_long_soft.h              pppm.cpp
action pair_lj_charmm_coul_long_soft.cpp  pppm.cpp
action pair_lj_charmm_coul_long_soft.h    pppm.cpp
action pair_lj_cut_coul_cut_soft.cpp
action pair_lj_cut_coul_cut_soft.h
action pair_lj_cut_coul_long_soft.cpp     pppm.cpp
action pair_lj_cut_coul_long_soft.h       pppm.cpp
action pair_lj_cut_soft.cpp
action pair_lj_cut_soft.h
action pair_lj_cut_tip4p_long_soft.cpp    pppm_tip4p.cpp
action pair_lj_cut_tip4p_long_soft.h      pppm_tip4p.cpp
action pair_morse_soft.cpp
action pair_morse_soft.h
action pair_tip4p_long_soft.cpp           pppm_tip4p.cpp
action pair_tip4p_long_soft.h             pppm_tip4p.cpp

