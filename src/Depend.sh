# Depend.sh = Install/unInstall files due to package dependencies
# this script is invoked after any package is installed/uninstalled

# parent package = has files that files in another package derive from
# child package = has files that derive from files in another package
# all parent/child package dependencies should be listed below

# re-install child packages that depend on the parent,
#   but only if the child package is already installed
# this is necessary to insure the child package installs
#   only child files whose parent package files are now installed
# decisions on (un)installing individual child files are made by
#   the Install.sh script in the child package

# depend function: args = child-package 0/1
# checks if child-package is installed, if not just return
# if parent package is being installed, reinstall the child
# if parent package is being uninstalled, uninstall the child, reinstall it

depend () {
  cd $1
  installed=0
  for file in *.cpp *.h; do
    if (test -e ../$file) then
      installed=1
    fi
  done
  cd ..
  if (test $installed = 0) then
    return
  fi

  if (test $2 = 1) then
    echo "  reinstalling package $1"
    cd $1; /bin/sh Install.sh 1; cd ..
  elif (test $2 = 0) then
    echo "  un/reinstalling package $1"
    cd $1; /bin/sh Install.sh 0; /bin/sh Install.sh 1; cd ..
  fi
}

# add one if statement per parent package
# add one depend() call per child package that depends on that parent
# list of child packages:
#   GPU, OPT, USER-CUDA, USER-MISC, USER-OMP

if (test $1 = "ASPHERE") then
  depend GPU $2
  depend USER-OMP $2
fi

if (test $1 = "CLASS2") then
  depend GPU $2
  depend USER-CUDA $2
  depend USER-OMP $2
fi

if (test $1 = "COLLOID") then
  depend GPU $2
  depend USER-OMP $2
fi

if (test $1 = "GRANULAR") then
  depend USER-CUDA $2
  depend USER-OMP $2
fi

if (test $1 = "KSPACE") then
  depend GPU $2
  depend OPT $2
  depend USER-CUDA $2
  depend USER-OMP $2
fi

if (test $1 = "MANYBODY") then
  depend GPU $2
  depend OPT $2
  depend USER-CUDA $2
  depend USER-MISC $2
  depend USER-OMP $2
fi

if (test $1 = "PERI") then
  depend USER-OMP $2
fi

if (test $1 = "RIGID") then
  depend USER-OMP $2
fi

if (test $1 = "USER-CG-CMM") then
  depend GPU $2
  depend USER-CUDA $2
  depend USER-OMP $2
fi

if (test $1 = "USER-MISC") then
  depend GPU $2
  depend USER-OMP $2
fi

