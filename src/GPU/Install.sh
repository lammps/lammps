# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude GPU info
# do not install child files if parent does not exist

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
include ..\/..\/lib\/gpu\/Makefile.lammps\
' ../Makefile.package.settings
  fi

  if (test -e ../pair_yukawa.cpp) then
    cp pair_yukawa_gpu.cpp ..
    cp pair_yukawa_gpu.h ..
  fi
  
  if (test -e ../pair_table.cpp) then
    cp pair_table_gpu.cpp ..
    cp pair_table_gpu.h ..
  fi

  if (test -e ../pair_buck.cpp) then
    cp pair_buck_gpu.cpp ..
    cp pair_buck_gpu.h ..
  fi

  if (test -e ../pair_buck_coul_cut.cpp) then
    cp pair_buck_coul_cut_gpu.cpp ..
    cp pair_buck_coul_cut_gpu.h ..
  fi

  if (test -e ../pair_buck_coul_long.cpp) then
    cp pair_buck_coul_long_gpu.cpp ..
    cp pair_buck_coul_long_gpu.h ..
  fi
  
  if (test -e ../pair_eam.cpp) then
    cp pair_eam_gpu.cpp ..
    cp pair_eam_gpu.h ..
    cp pair_eam_alloy_gpu.cpp ..
    cp pair_eam_alloy_gpu.h ..
    cp pair_eam_fs_gpu.cpp ..
    cp pair_eam_fs_gpu.h ..
  fi
  
  if (test -e ../pair_gayberne.cpp) then
    cp pair_gayberne_gpu.cpp ..
    cp pair_gayberne_gpu.h ..
    cp pair_resquared_gpu.cpp ..
    cp pair_resquared_gpu.h ..
  fi
  
  if (test -e ../pair_lj_cut_coul_long.cpp) then
    cp pair_lj_cut_coul_long_gpu.cpp ..
    cp pair_lj_cut_coul_long_gpu.h ..
  fi

  if (test -e ../pair_lj_class2.cpp) then
    cp pair_lj_class2_gpu.cpp ..
    cp pair_lj_class2_gpu.h ..
  fi

  if (test -e ../pair_lj_class2_coul_long.cpp) then
    cp pair_lj_class2_coul_long_gpu.cpp ..
    cp pair_lj_class2_coul_long_gpu.h ..
  fi

  if (test -e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_lj_charmm_coul_long_gpu.cpp ..
    cp pair_lj_charmm_coul_long_gpu.h ..
  fi

  if (test -e ../pair_coul_long.cpp) then
    cp pair_coul_long_gpu.cpp ..
    cp pair_coul_long_gpu.h ..
  fi

  if (test -e ../pair_lj_sdk.cpp) then
    cp pair_lj_sdk_gpu.cpp ..
    cp pair_lj_sdk_gpu.h ..
  fi

  if (test -e ../pair_lj_sdk_coul_long.cpp) then
    cp pair_lj_sdk_coul_long_gpu.cpp ..
    cp pair_lj_sdk_coul_long_gpu.h ..
  fi

  if (test -e ../pppm.cpp) then
    cp pppm_gpu.cpp ..
    cp pppm_gpu.h ..
  fi

  cp pair_lj_cut_gpu.cpp ..
  cp pair_morse_gpu.cpp ..
  cp pair_lj96_cut_gpu.cpp ..
  cp pair_lj_expand_gpu.cpp ..
  cp pair_lj_cut_coul_cut_gpu.cpp ..

  cp fix_gpu.cpp ..

  cp pair_lj_cut_gpu.h ..
  cp pair_morse_gpu.h ..
  cp pair_lj96_cut_gpu.h ..
  cp pair_lj_expand_gpu.h ..
  cp pair_lj_cut_coul_cut_gpu.h ..
  
  cp fix_gpu.h ..
  cp gpu_extra.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu[^ \t]* //' ../Makefile.package
  fi
  
  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*gpu.*$/d' ../Makefile.package.settings
  fi

  rm -f ../pair_buck_coul_cut_gpu.cpp
  rm -f ../pair_buck_coul_long_gpu.cpp
  rm -f ../pair_buck_gpu.cpp
  rm -f ../pair_coul_long_gpu.cpp
  rm -f ../pair_eam_alloy_gpu.cpp
  rm -f ../pair_eam_fs_gpu.cpp
  rm -f ../pair_eam_gpu.cpp
  rm -f ../pair_gayberne_gpu.cpp
  rm -f ../pair_lj96_cut_gpu.cpp
  rm -f ../pair_lj_charmm_coul_long_gpu.cpp
  rm -f ../pair_lj_class2_coul_long_gpu.cpp
  rm -f ../pair_lj_class2_gpu.cpp
  rm -f ../pair_lj_cut_coul_cut_gpu.cpp
  rm -f ../pair_lj_cut_coul_long_gpu.cpp
  rm -f ../pair_lj_cut_gpu.cpp
  rm -f ../pair_lj_expand_gpu.cpp
  rm -f ../pair_lj_sdk_coul_long_gpu.cpp
  rm -f ../pair_lj_sdk_gpu.cpp
  rm -f ../pair_morse_gpu.cpp
  rm -f ../pair_resquared_gpu.cpp
  rm -f ../pair_table_gpu.cpp
  rm -f ../pair_yukawa_gpu.cpp
  rm -f ../pppm_gpu.cpp

  rm -f ../fix_gpu.cpp

  rm -f ../pair_buck_coul_cut_gpu.h
  rm -f ../pair_buck_coul_long_gpu.h
  rm -f ../pair_buck_gpu.h
  rm -f ../pair_coul_long_gpu.h
  rm -f ../pair_eam_alloy_gpu.h
  rm -f ../pair_eam_fs_gpu.h
  rm -f ../pair_eam_gpu.h
  rm -f ../pair_gayberne_gpu.h
  rm -f ../pair_lj96_cut_gpu.h
  rm -f ../pair_lj_charmm_coul_long_gpu.h
  rm -f ../pair_lj_class2_coul_long_gpu.h
  rm -f ../pair_lj_class2_gpu.h
  rm -f ../pair_lj_cut_coul_cut_gpu.h
  rm -f ../pair_lj_cut_coul_long_gpu.h
  rm -f ../pair_lj_cut_gpu.h
  rm -f ../pair_lj_expand_gpu.h
  rm -f ../pair_lj_sdk_coul_long_gpu.h
  rm -f ../pair_lj_sdk_gpu.h
  rm -f ../pair_morse_gpu.h
  rm -f ../pair_resquared_gpu.h
  rm -f ../pair_table_gpu.h
  rm -f ../pair_yukawa_gpu.h
  rm -f ../pppm_gpu.h

  rm -f ../fix_gpu.h
  rm -f ../gpu_extra.h
  
fi

