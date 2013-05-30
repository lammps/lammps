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
  elif (test ! -z "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# list of files with optional dependencies

action angle_cosine_shift.cpp
action angle_cosine_shift.h
action angle_cosine_shift_exp.cpp
action angle_cosine_shift_exp.h
action angle_dipole.cpp
action angle_dipole.h
action angle_fourier.cpp
action angle_fourier.h
action angle_fourier_simple.cpp
action angle_fourier_simple.h
action angle_quartic.cpp
action angle_quartic.h
action bond_harmonic_shift.cpp
action bond_harmonic_shift.h
action bond_harmonic_shift_cut.cpp
action bond_harmonic_shift_cut.h
action compute_ackland_atom.cpp
action compute_ackland_atom.h
action compute_temp_rotate.cpp
action compute_temp_rotate.h
action dihedral_cosine_shift_exp.cpp
action dihedral_cosine_shift_exp.h
action dihedral_fourier.cpp
action dihedral_fourier.h
action dihedral_nharmonic.cpp
action dihedral_nharmonic.h
action dihedral_quadratic.cpp
action dihedral_quadratic.h
action dihedral_table.cpp
action dihedral_table.h
action fix_addtorque.cpp
action fix_addtorque.h
action fix_imd.cpp
action fix_imd.h
action fix_smd.cpp
action fix_smd.h
action improper_cossq.cpp
action improper_cossq.h
action improper_fourier.cpp
action improper_fourier.h
action improper_ring.cpp
action improper_ring.h
action pair_cdeam.cpp pair_eam_alloy.cpp
action pair_cdeam.h pair_eam_alloy.cpp
action pair_coul_diel.cpp
action pair_coul_diel.h
action pair_dipole_sf.cpp
action pair_dipole_sf.h
action pair_edip.cpp
action pair_edip.h
action pair_gauss_cut.cpp
action pair_gauss_cut.h
action pair_lj_sf.cpp
action pair_lj_sf.h
action pair_meam_spline.cpp
action pair_meam_spline.h
action pair_meam_sw_spline.cpp
action pair_meam_sw_spline.h
action pair_tersoff_table.cpp
action pair_tersoff_table.h
