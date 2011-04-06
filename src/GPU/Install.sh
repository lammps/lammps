# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude GPU library
# do not copy gayberne files if non-GPU version does not exist
# do not copy charmm files if non-GPU version does not exist

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
    sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L../../lib/gpu |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lgpu |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(gpu_SYSPATH) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(gpu_SYSLIB) |' ../Makefile.package
  fi
  
  if (test -e ../pppm.cpp) then
    cp pppm_gpu.cpp ..
    cp pppm_gpu_single.cpp ..
    cp pppm_gpu_double.cpp ..
    cp pppm_gpu.h ..
    cp pppm_gpu_single.h ..
    cp pppm_gpu_double.h ..
  fi
  
  if (test -e ../pair_gayberne.cpp) then
    cp pair_gayberne_gpu.cpp ..
    cp pair_gayberne_gpu.h ..
  fi
  
  if (test -e ../pair_lj_cut_coul_long.cpp) then
    cp pair_lj_cut_coul_long_gpu.cpp ..
    cp pair_lj_cut_coul_long_gpu.h ..
  fi

  if (test -e ../pair_lj_charmm_coul_long.cpp) then
    cp pair_lj_charmm_coul_long_gpu.cpp ..
    cp pair_lj_charmm_coul_long_gpu.h ..
    cp pair_lj_charmm_coul_long_gpu2.cpp ..
    cp pair_lj_charmm_coul_long_gpu2.h ..
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

  cp pair_lj_cut_gpu.cpp ..
  cp pair_morse_gpu.cpp ..
  cp pair_lj96_cut_gpu.cpp ..
  cp pair_lj_expand_gpu.cpp ..
  cp pair_lj_cut_coul_cut_gpu.cpp ..
  cp pair_lj_cut_tgpu.cpp ..
  cp pair_lj_cut_gpu.h ..
  cp pair_morse_gpu.h ..
  cp pair_lj96_cut_gpu.h ..
  cp pair_lj_expand_gpu.h ..
  cp pair_lj_cut_coul_cut_gpu.h ..
  cp pair_lj_cut_tgpu.h ..
  
  cp fix_gpu.cpp ..
  cp fix_gpu.h ..
  cp gpu_extra.h ..

  cp pair_omp_gpu.cpp ..
  cp pair_omp_gpu.h ..
  
elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*gpu //' ../Makefile.package
    sed -i -e 's/[^ \t]*gpu_[^ \t]*) //' ../Makefile.package
  fi
  
  rm ../pppm_gpu.cpp
  rm ../pppm_gpu_single.cpp
  rm ../pppm_gpu_double.cpp
  rm ../pair_gayberne_gpu.cpp
  rm ../pair_lj_cut_gpu.cpp
  rm ../pair_morse_gpu.cpp
  rm ../pair_lj96_cut_gpu.cpp
  rm ../pair_lj_expand_gpu.cpp
  rm ../pair_lj_cut_coul_cut_gpu.cpp
  rm ../pair_lj_cut_coul_long_gpu.cpp
  rm ../pair_lj_charmm_coul_long_gpu.cpp
  rm ../pair_lj_charmm_coul_long_gpu2.cpp
  rm ../pair_lj_cut_tgpu.cpp
  rm ../pair_cg_cmm_gpu.cpp
  rm ../pair_cg_cmm_coul_long_gpu.cpp
  rm ../pair_cg_cmm_coul_msm.cpp
  rm ../pair_cg_cmm_coul_msm_gpu.cpp
  rm ../fix_gpu.cpp
  rm ../pair_omp_gpu.cpp

  rm ../pppm_gpu.h
  rm ../pppm_gpu_single.cpp
  rm ../pppm_gpu_double.h
  rm ../pair_gayberne_gpu.h
  rm ../pair_lj_cut_gpu.h
  rm ../pair_morse_gpu.h
  rm ../pair_lj96_cut_gpu.h
  rm ../pair_lj_expand_gpu.h
  rm ../pair_lj_cut_coul_cut_gpu.h
  rm ../pair_lj_cut_coul_long_gpu.h
  rm ../pair_lj_charmm_coul_long_gpu.h
  rm ../pair_lj_charmm_coul_long_gpu2.h
  rm ../pair_lj_cut_tgpu.cpp
  rm ../pair_cg_cmm_gpu.h
  rm ../pair_cg_cmm_coul_long_gpu.h
  rm ../pair_cg_cmm_coul_msm.h
  rm ../pair_cg_cmm_coul_msm_gpu.h
  rm ../fix_gpu.h
  rm ../gpu_extra.h
  rm ../pair_omp_gpu.h
  
fi

