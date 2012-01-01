# Depend.sh = Install/unInstall files for dependent packages
# only Install/unInstall if dependent package is already installed
# all packages with dependencies should be listed here
# install dependent child files when parent files installed
# uninstall dependent child files when parent files uninstalled
# decisions on individual files are made by package Install.sh scripts

if (test $1 = 1) then

  if (test -e pair_lj_cut_opt.h) then
    cd OPT; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_gpu.h) then
    cd GPU; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e cg_cmm_params.h) then
    cd USER-CG-CMM; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_cuda.h) then
    cd USER-CUDA; /bin/sh Install.sh 1; cd ..
  fi

elif (test $1 = 0) then

  if (test -e pair_lj_cut_opt.h) then
    cd OPT; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_gpu.h) then
    cd GPU; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e cg_cmm_params.h) then
    cd USER-CG-CMM; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
  if (test -e pair_lj_cut_cuda.h) then
    cd USER-CUDA; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi

fi
