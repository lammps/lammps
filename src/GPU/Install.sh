# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

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

action fix_gpu.cpp
action fix_gpu.h
action gpu_extra.h
action pair_beck_gpu.cpp
action pair_beck_gpu.h 
action pair_born_coul_long_gpu.cpp pair_born_coul_long.cpp
action pair_born_coul_long_gpu.h pair_born_coul_long.cpp
action pair_born_coul_wolf_gpu.cpp
action pair_born_coul_wolf_gpu.h
action pair_born_gpu.cpp
action pair_born_gpu.h
action pair_buck_coul_cut_gpu.cpp pair_buck_coul_cut.cpp
action pair_buck_coul_cut_gpu.h pair_buck_coul_cut.cpp
action pair_buck_coul_long_gpu.cpp pair_buck_coul_long.cpp
action pair_buck_coul_long_gpu.h pair_buck_coul_long.cpp
action pair_buck_gpu.cpp pair_buck.cpp
action pair_buck_gpu.h pair_buck.cpp
action pair_colloid_gpu.cpp pair_colloid.cpp
action pair_colloid_gpu.h pair_colloid.cpp
action pair_coul_dsf_gpu.cpp
action pair_coul_dsf_gpu.h
action pair_coul_long_gpu.cpp pair_coul_long.cpp
action pair_coul_long_gpu.h pair_coul_long.cpp
action pair_lj_cut_dipole_cut_gpu.cpp pair_lj_cut_dipole_cut.cpp
action pair_lj_cut_dipole_cut_gpu.h pair_lj_cut_dipole_cut.cpp
action pair_lj_sf_dipole_sf_gpu.cpp pair_lj_sf_dipole_sf.cpp
action pair_lj_sf_dipole_sf_gpu.h pair_lj_sf_dipole_sf.cpp
action pair_eam_alloy_gpu.cpp pair_eam.cpp
action pair_eam_alloy_gpu.h pair_eam.cpp
action pair_eam_fs_gpu.cpp pair_eam.cpp 
action pair_eam_fs_gpu.h pair_eam.cpp
action pair_eam_gpu.cpp pair_eam.cpp 
action pair_eam_gpu.h pair_eam.cpp
action pair_gauss_gpu.cpp
action pair_gauss_gpu.h
action pair_gayberne_gpu.cpp pair_gayberne.cpp
action pair_gayberne_gpu.h pair_gayberne.cpp
action pair_lj96_cut_gpu.cpp
action pair_lj96_cut_gpu.h
action pair_lj_charmm_coul_long_gpu.cpp pair_lj_charmm_coul_long.cpp
action pair_lj_charmm_coul_long_gpu.h pair_lj_charmm_coul_long.cpp
action pair_lj_class2_coul_long_gpu.cpp pair_lj_class2_coul_long.cpp
action pair_lj_class2_coul_long_gpu.h pair_lj_class2_coul_long.cpp
action pair_lj_class2_gpu.cpp pair_lj_class2.cpp
action pair_lj_class2_gpu.h pair_lj_class2.cpp
action pair_lj_cut_coul_cut_gpu.cpp
action pair_lj_cut_coul_cut_gpu.h
action pair_lj_cut_coul_debye_gpu.cpp
action pair_lj_cut_coul_debye_gpu.h
action pair_lj_cut_coul_dsf_gpu.cpp
action pair_lj_cut_coul_dsf_gpu.h
action pair_lj_cut_coul_long_gpu.cpp pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_long_gpu.h pair_lj_cut_coul_long.cpp
action pair_lj_cut_coul_msm_gpu.cpp
action pair_lj_cut_coul_msm_gpu.h
action pair_lj_cut_gpu.cpp
action pair_lj_cut_gpu.h
action pair_lj_expand_gpu.cpp
action pair_lj_expand_gpu.h
action pair_lj_gromacs_gpu.cpp
action pair_lj_gromacs_gpu.h
action pair_lj_sdk_coul_long_gpu.cpp pair_lj_sdk_coul_long.cpp
action pair_lj_sdk_coul_long_gpu.h pair_lj_sdk_coul_long.cpp
action pair_lj_sdk_gpu.cpp pair_lj_sdk.cpp
action pair_lj_sdk_gpu.h pair_lj_sdk.cpp
action pair_mie_cut_gpu.cpp
action pair_mie_cut_gpu.h 
action pair_morse_gpu.cpp
action pair_morse_gpu.h
action pair_resquared_gpu.cpp pair_resquared.cpp
action pair_resquared_gpu.h pair_resquared.cpp
action pair_soft_gpu.cpp
action pair_soft_gpu.h
action pair_sw_gpu.cpp pair_sw.cpp
action pair_sw_gpu.h pair_sw.h
action pair_table_gpu.cpp pair_table.cpp
action pair_table_gpu.h pair_table.cpp
action pair_yukawa_colloid_gpu.cpp pair_yukawa_colloid.cpp
action pair_yukawa_colloid_gpu.h pair_yukawa_colloid.cpp
action pair_yukawa_gpu.cpp pair_yukawa.cpp
action pair_yukawa_gpu.h pair_yukawa.cpp
action pppm_gpu.cpp pppm.cpp
action pppm_gpu.h pppm.cpp

# edit 2 Makefile.package files to include/exclude GPU info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(gpu_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*gpu.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/gpu\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu[^ \t]* //' ../Makefile.package
  fi
  
  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*gpu.*$/d' ../Makefile.package.settings
  fi

fi
