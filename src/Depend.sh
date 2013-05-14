# Depend.sh = Install/unInstall files for dependent packages
# all packages which contain one or more files that depend on
#   other packages should be listed here, in both the 1 and 0 clauses
# this script is invoked after any parent package is installed/uninstalled
# this script re-installs child packages that depend on the parent,
#   but only if the child package is already installed
# this is necessary to insure the child package installs
#   only child files whose parent package files are now installed
# decisions on installing individual child files are made by
#   the Install.sh script in the child package

if (test $1 = 1) then

  if (test -e pair_lj_cut_gpu.h) then
    cd GPU; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_opt.h) then
    cd OPT; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e cg_cmm_params.h) then
    cd USER-CG-CMM; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_cuda.h) then
    cd USER-CUDA; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e fix_imd.h) then
    cd USER-MISC; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e thr_omp.h) then
    cd USER-OMP; /bin/sh Install.sh 1; cd ..
  fi

elif (test $1 = 0) then

  if (test -e pair_lj_cut_gpu.h) then
    cd GPU; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_opt.h) then
    cd OPT; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e cg_cmm_params.h) then
    cd USER-CG-CMM; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_cuda.h) then
    cd USER-CUDA; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e fix_imd.h) then
    cd USER-MISC; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e thr_omp.h) then
    cd USER-OMP; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi

fi
