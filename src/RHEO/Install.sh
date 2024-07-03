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

# package files without dependencies
action atom_vec_rheo_thermal.h
action atom_vec_rheo_thermal.cpp
action atom_vec_rheo.h
action atom_vec_rheo.cpp
action compute_rheo_grad.h
action compute_rheo_grad.cpp
action compute_rheo_interface.h
action compute_rheo_interface.cpp
action compute_rheo_kernel.h
action compute_rheo_kernel.cpp
action compute_rheo_rho_sum.h
action compute_rheo_rho_sum.cpp
action compute_rheo_surface.h
action compute_rheo_surface.cpp
action compute_rheo_vshift.h
action compute_rheo_vshift.cpp
action fix_rheo_oxidation.h
action fix_rheo_oxidation.cpp
action fix_rheo_pressure.h
action fix_rheo_pressure.cpp
action fix_rheo_viscosity.h
action fix_rheo_viscosity.cpp
action fix_rheo.h
action fix_rheo.cpp
action pair_rheo.h
action pair_rheo.cpp
action pair_rheo_solid.h
action pair_rheo_solid.cpp

# package files with dependencies
action bond_rheo_shell.h bond_bpm.h
action bond_rheo_shell.cpp bond_bpm.h
action compute_rheo_property_atom.h fix_update_special_bonds.h
action compute_rheo_property_atom.cpp fix_update_special_bonds.h
action fix_rheo_thermal.h fix_update_special_bonds.h
action fix_rheo_thermal.cpp fix_update_special_bonds.h

# Warn that some styles in RHEO have base classes in BPM

if (test $1 = 1) then
  if (test ! -e ../bond_bpm.cpp) then
    echo "Must install BPM package to use all features of RHEO package"
  fi
fi

# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*rheo[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(rheo_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(rheo_SYSLIB) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*rheo.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/..\/lib\/rheo\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*rheo[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*rheo.*$/d' ../Makefile.package.settings
  fi

fi
