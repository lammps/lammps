# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL

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

action pair_eam_alloy_opt.cpp pair_eam_alloy.cpp
action pair_eam_alloy_opt.h pair_eam_alloy.cpp
action pair_eam_fs_opt.cpp pair_eam_fs.cpp
action pair_eam_fs_opt.h pair_eam_fs.cpp
action pair_eam_opt.cpp pair_eam.cpp
action pair_eam_opt.h pair_eam.cpp
action pair_lj_charmm_coul_long_opt.cpp pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_long_opt.h pair_lj_charmm_coul_long.cpp
action pair_lj_cut_coul_long_opt.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_opt.h pair_lj_cut_coul_long.cpp
action pair_lj_cut_opt.cpp
action pair_lj_cut_opt.h
action pair_lj_cut_tip4p_long_opt.cpp pair_lj_cut_tip4p_long.cpp
action pair_lj_cut_tip4p_long_opt.h pair_lj_cut_tip4p_long.cpp
action pair_lj_long_coul_long_opt.cpp pair_lj_long_coul_long.cpp
action pair_lj_long_coul_long_opt.h pair_lj_long_coul_long.cpp
action pair_morse_opt.cpp
action pair_morse_opt.h
action pair_ufm_opt.cpp pair_ufm.cpp
action pair_ufm_opt.h pair_ufm.h
action pair_ilp_graphene_hbn_opt.cpp pair_ilp_graphene_hbn.cpp
action pair_ilp_graphene_hbn_opt.h pair_ilp_graphene_hbn.h
action pair_ilp_tmd_opt.cpp pair_ilp_tmd.cpp
action pair_ilp_tmd_opt.h pair_ilp_tmd.h
action pair_saip_metal_opt.cpp pair_saip_metal.cpp
action pair_saip_metal_opt.h pair_saip_metal.h
