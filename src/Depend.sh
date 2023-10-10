# Depend.sh = Install/unInstall files due to package dependencies
# this script is invoked after any package is installed/uninstalled

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# all parent/child package dependencies should be listed below
# parent package = has files that files in another package derive from
# child package = has files that derive from files in another package

# update child packages that depend on the parent,
#   but only if the child package is already installed
# this is necessary to ensure the child package installs
#   only child files whose parent package files are now installed
# decisions on (un)installing individual child files are made by
#   the Install.sh script in the child package

# depend function: arg = child-package
# checks if child-package is installed, if not just return
# otherwise invoke update of child package via its Install.sh

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

  echo "  updating package $1"
  if (test -e $1/Install.sh) then
    cd $1; /bin/sh Install.sh 2; cd ..
  else
    cd $1; /bin/sh ../Install.sh 2; cd ..
  fi
}

# add one if statement per parent package
# add one depend() call per child package that depends on that parent

if (test $1 = "AMOEBA") then
  depend GPU
fi

if (test $1 = "ASPHERE") then
  depend GPU
  depend OPENMP
  depend CG-DNA
  depend INTEL
fi

if (test $1 = "CLASS2") then
  depend GPU
  depend KOKKOS
  depend OPENMP
fi

if (test $1 = "COLLOID") then
  depend GPU
  depend OPENMP
fi

if (test $1 = "DIELECTRIC") then
  depend OPENMP
fi

if (test $1 = "DIPOLE") then
  depend KOKKOS
  depend OPENMP
fi

if (test $1 = "DPD-BASIC") then
  depend GPU
  depend KOKKOS
  depend OPENMP
  depend INTEL
fi

if (test $1 = "ELECTRODE") then
  depend INTEL
fi

if (test $1 = "EXTRA-COMPUTE") then
  depend KOKKOS
fi

if (test $1 = "EXTRA-MOLECULE") then
  depend GPU
  depend OPENMP
fi

if (test $1 = "EXTRA-PAIR") then
  depend GPU
  depend OPENMP
fi

if (test $1 = "GRANULAR") then
  depend KOKKOS
  depend OPENMP
fi

if (test $1 = "INTERLAYER") then
  depend OPT
fi

if (test $1 = "KSPACE") then
  depend CG-SPICA
  depend CORESHELL
  depend DIELECTRIC
  depend GPU
  depend KOKKOS
  depend OPT
  depend OPENMP
  depend INTEL
  depend PHONON
  depend FEP
fi

if (test $1 = "LEPTON") then
  depend OPENMP
fi

if (test $1 = "MANYBODY") then
  depend ATC
  depend GPU
  depend INTEL
  depend KOKKOS
  depend OPT
  depend QEQ
  depend OPENMP
fi

if (test $1 = "MC") then
  depend MISC
fi

if (test $1 = "MEAM") then
  depend KOKKOS
fi

if (test $1 = "MOLECULE") then
  depend EXTRA-MOLECULE
  depend GPU
  depend KOKKOS
  depend FEP
  depend MESONT
  depend OPENMP
  depend INTEL
fi

if (test $1 = "PERI") then
  depend OPENMP
fi

if (test $1 = "PYTHON") then
  depend ML-IAP
  depend KOKKOS
fi

if (test $1 = "PHONON") then
  depend KOKKOS
fi

if (test $1 = "RIGID") then
  depend KOKKOS
  depend OPENMP
  depend DPD-SMOOTH
fi

if (test $1 = "ML-IAP") then
  depend KOKKOS
fi

if (test $1 = "ML-PACE") then
  depend KOKKOS
fi

if (test $1 = "ML-SNAP") then
  depend ML-IAP
  depend KOKKOS
  depend INTEL
fi

if (test $1 = "CG-SPICA") then
  depend GPU
  depend KOKKOS
  depend OPENMP
fi

if (test $1 = "DPD-REACT") then
  depend KOKKOS
fi

if (test $1 = "DRUDE") then
  depend OPENMP
fi

if (test $1 = "FEP") then
  depend OPENMP
fi

if (test $1 = "REAXFF") then
  depend KOKKOS
  depend OPENMP
fi

if (test $1 = "SPIN") then
  depend KOKKOS
fi
