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

# list of files with optional dependcies

action atom_vec_dielectric.cpp
action atom_vec_dielectric.h
action fix_polarize_bem_gmres.cpp
action fix_polarize_bem_gmres.h
action fix_polarize_bem_icc.cpp
action fix_polarize_bem_icc.h
action fix_polarize_functional.cpp
action fix_polarize_functional.h
action pair_lj_cut_coul_msm_dielectric.cpp
action pair_lj_cut_coul_msm_dielectric.h
action pair_lj_cut_coul_long_dielectric.cpp
action pair_lj_cut_coul_long_dielectric.h
action pair_lj_long_coul_long_dielectric.cpp pair_lj_long_coul_long.cpp
action pair_lj_long_coul_long_dielectric.h pair_lj_long_coul_long.cpp
action pair_lj_cut_coul_cut_dielectric.cpp
action pair_lj_cut_coul_cut_dielectric.h
action pair_coul_long_dielectric.cpp
action pair_coul_long_dielectric.h
action pair_coul_cut_dielectric.cpp
action pair_coul_cut_dielectric.h
action pppm_dielectric.cpp
action pppm_dielectric.h
action pppm_disp_dielectric.cpp pppm_disp.cpp
action pppm_disp_dielectric.h pppm_disp.h
action msm_dielectric.cpp
action msm_dielectric.h
