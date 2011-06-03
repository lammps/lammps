# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude GPU info
# do not install child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
    sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(gpu_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
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

  if (test -e ../pair_cg_cmm.cpp) then
    cp pair_cg_cmm_gpu.cpp ..
    cp pair_cg_cmm_gpu.h ..
  fi

  if (test -e ../pair_cg_cmm_coul_long.cpp) then
    cp pair_cg_cmm_coul_long_gpu.cpp ..
    cp pair_cg_cmm_coul_long_gpu.h ..
    cp pair_cg_cmm_coul_msm.cpp ..
    cp pair_cg_cmm_coul_msm.h ..
    cp pair_cg_cmm_coul_msm_gpu.cpp ..
    cp pair_cg_cmm_coul_msm_gpu.h ..
  fi

  if (test -e ../pppm.cpp) then
    cp pppm_gpu.cpp ..
    cp pppm_gpu_single.cpp ..
    cp pppm_gpu_double.cpp ..
    cp pppm_gpu.h ..
    cp pppm_gpu_single.h ..
    cp pppm_gpu_double.h ..
  fi

  cp pair_lj_cut_gpu.cpp ..
  cp pair_morse_gpu.cpp ..
  cp pair_lj96_cut_gpu.cpp ..
  cp pair_lj_expand_gpu.cpp ..
  cp pair_lj_cut_coul_cut_gpu.cpp ..
  cp pair_lj_cut_tgpu.cpp ..

  cp fix_gpu.cpp ..

  cp pair_lj_cut_gpu.h ..
  cp pair_morse_gpu.h ..
  cp pair_lj96_cut_gpu.h ..
  cp pair_lj_expand_gpu.h ..
  cp pair_lj_cut_coul_cut_gpu.h ..
  cp pair_lj_cut_tgpu.h ..
  
  cp fix_gpu.h ..
  cp gpu_extra.h ..

  cp pair_omp_gpu.cpp ..
  cp pair_omp_gpu.h ..
  
elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
    sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
  fi
  
  rm -f ../pppm_gpu.cpp
  rm -f ../pppm_gpu_single.cpp
  rm -f ../pppm_gpu_double.cpp
  rm -f ../pair_gayberne_gpu.cpp
  rm -f ../pair_resquared_gpu.cpp
  rm -f ../pair_lj_cut_gpu.cpp
  rm -f ../pair_morse_gpu.cpp
  rm -f ../pair_lj96_cut_gpu.cpp
  rm -f ../pair_lj_expand_gpu.cpp
  rm -f ../pair_lj_cut_coul_cut_gpu.cpp
  rm -f ../pair_lj_cut_coul_long_gpu.cpp
  rm -f ../pair_lj_class2_gpu.cpp
  rm -f ../pair_lj_class2_coul_long_gpu.cpp
  rm -f ../pair_lj_charmm_coul_long_gpu.cpp
  rm -f ../pair_lj_cut_tgpu.cpp
  rm -f ../pair_cg_cmm_gpu.cpp
  rm -f ../pair_cg_cmm_coul_long_gpu.cpp
  rm -f ../pair_cg_cmm_coul_msm.cpp
  rm -f ../pair_cg_cmm_coul_msm_gpu.cpp

  rm -f ../fix_gpu.cpp
  rm -f ../pair_omp_gpu.cpp

  rm -f ../pppm_gpu.h
  rm -f ../pppm_gpu_single.h
  rm -f ../pppm_gpu_double.h
  rm -f ../pair_gayberne_gpu.h
  rm -f ../pair_resquared_gpu.h
  rm -f ../pair_lj_cut_gpu.h
  rm -f ../pair_morse_gpu.h
  rm -f ../pair_lj96_cut_gpu.h
  rm -f ../pair_lj_expand_gpu.h
  rm -f ../pair_lj_cut_coul_cut_gpu.h
  rm -f ../pair_lj_cut_coul_long_gpu.h
  rm -f ../pair_lj_class2_gpu.h
  rm -f ../pair_lj_class2_coul_long_gpu.h
  rm -f ../pair_lj_charmm_coul_long_gpu.h
  rm -f ../pair_lj_cut_tgpu.cpp
  rm -f ../pair_cg_cmm_gpu.h
  rm -f ../pair_cg_cmm_coul_long_gpu.h
  rm -f ../pair_cg_cmm_coul_msm.h
  rm -f ../pair_cg_cmm_coul_msm_gpu.h

  rm -f ../fix_gpu.h
  rm -f ../gpu_extra.h
  rm -f ../pair_omp_gpu.h
  
fi

