# Install/unInstall package classes in LAMMPS

# bond.h, angle.h, dihedal.h, improper.h must always be in src
# bond_hybrid.h must always be in src

if ($1 == 1) then

  cp style_molecule.h ..

  cp angle.cpp ..
  cp angle_charmm.cpp ..
  cp angle_cosine.cpp ..
  cp angle_cosine_delta.cpp ..
  cp angle_cosine_squared.cpp ..
  cp angle_harmonic.cpp ..
  cp angle_hybrid.cpp ..
  cp angle_table.cpp ..
  cp atom_vec_angle.cpp ..
  cp atom_vec_bond.cpp ..
  cp atom_vec_full.cpp ..
  cp atom_vec_molecular.cpp ..
  cp bond.cpp ..
  cp bond_fene.cpp ..
  cp bond_fene_expand.cpp ..
  cp bond_harmonic.cpp ..
  cp bond_hybrid.cpp ..
  cp bond_morse.cpp ..
  cp bond_nonlinear.cpp ..
  cp bond_quartic.cpp ..
  cp bond_table.cpp ..
  cp dihedral.cpp ..
  cp dihedral_charmm.cpp ..
  cp dihedral_harmonic.cpp ..
  cp dihedral_helix.cpp ..
  cp dihedral_hybrid.cpp ..
  cp dihedral_multi_harmonic.cpp ..
  cp dihedral_opls.cpp ..
  cp dump_bond.cpp ..
  cp fix_bond_break.cpp ..
  cp fix_bond_create.cpp ..
  cp fix_bond_swap.cpp ..
  cp improper.cpp ..
  cp improper_cvff.cpp ..
  cp improper_harmonic.cpp ..
  cp improper_hybrid.cpp ..
  cp pair_lj_charmm_coul_charmm.cpp ..
  cp pair_lj_charmm_coul_charmm_implicit.cpp ..

#  cp angle.h ..
  cp angle_charmm.h ..
  cp angle_cosine.h ..
  cp angle_cosine_delta.h ..
  cp angle_cosine_squared.h ..
  cp angle_harmonic.h ..
  cp angle_hybrid.h ..
  cp angle_table.h ..
  cp atom_vec_angle.h ..
  cp atom_vec_bond.h ..
  cp atom_vec_full.h ..
  cp atom_vec_molecular.h ..
#  cp bond.h ..
  cp bond_fene.h ..
  cp bond_fene_expand.h ..
  cp bond_harmonic.h ..
  cp bond_morse.h ..
#  cp bond_hybrid.h ..
  cp bond_nonlinear.h ..
  cp bond_quartic.h ..
  cp bond_table.h ..
#  cp dihedral.h ..
  cp dihedral_charmm.h ..
  cp dihedral_harmonic.h ..
  cp dihedral_helix.h ..
  cp dihedral_hybrid.h ..
  cp dihedral_multi_harmonic.h ..
  cp dihedral_opls.h ..
  cp dump_bond.h ..
  cp fix_bond_break.h ..
  cp fix_bond_create.h ..
  cp fix_bond_swap.h ..
#  cp improper.h ..
  cp improper_cvff.h ..
  cp improper_harmonic.h ..
  cp improper_hybrid.h ..
  cp pair_lj_charmm_coul_charmm.h ..
  cp pair_lj_charmm_coul_charmm_implicit.h ..

else if ($1 == 0) then

  rm ../style_molecule.h
  touch ../style_molecule.h

  rm ../angle.cpp
  rm ../angle_charmm.cpp
  rm ../angle_cosine.cpp
  rm ../angle_cosine_delta.cpp
  rm ../angle_cosine_squared.cpp
  rm ../angle_harmonic.cpp
  rm ../angle_hybrid.cpp
  rm ../angle_table.cpp
  rm ../atom_vec_angle.cpp
  rm ../atom_vec_bond.cpp
  rm ../atom_vec_full.cpp
  rm ../atom_vec_molecular.cpp
  rm ../bond.cpp
  rm ../bond_fene.cpp
  rm ../bond_fene_expand.cpp
  rm ../bond_harmonic.cpp
  rm ../bond_hybrid.cpp
  rm ../bond_morse.cpp
  rm ../bond_nonlinear.cpp
  rm ../bond_quartic.cpp
  rm ../bond_table.cpp
  rm ../dihedral.cpp
  rm ../dihedral_charmm.cpp
  rm ../dihedral_harmonic.cpp
  rm ../dihedral_helix.cpp
  rm ../dihedral_hybrid.cpp
  rm ../dihedral_multi_harmonic.cpp
  rm ../dihedral_opls.cpp
  rm ../dump_bond.cpp
  rm ../fix_bond_break.cpp
  rm ../fix_bond_create.cpp
  rm ../fix_bond_swap.cpp
  rm ../improper.cpp
  rm ../improper_cvff.cpp
  rm ../improper_harmonic.cpp
  rm ../improper_hybrid.cpp
  rm ../pair_lj_charmm_coul_charmm.cpp
  rm ../pair_lj_charmm_coul_charmm_implicit.cpp

#  rm ../angle.h
  rm ../angle_charmm.h
  rm ../angle_cosine.h
  rm ../angle_cosine_delta.h
  rm ../angle_cosine_squared.h
  rm ../angle_harmonic.h
  rm ../angle_hybrid.h
  rm ../angle_table.h
  rm ../atom_vec_angle.h
  rm ../atom_vec_bond.h
  rm ../atom_vec_full.h
  rm ../atom_vec_molecular.h
#  rm ../bond.h
  rm ../bond_fene.h
  rm ../bond_fene_expand.h
  rm ../bond_harmonic.h
#  rm ../bond_hybrid.h
  rm ../bond_morse.h
  rm ../bond_nonlinear.h
  rm ../bond_quartic.h
  rm ../bond_table.h
#  rm ../dihedral.h
  rm ../dihedral_charmm.h
  rm ../dihedral_harmonic.h
  rm ../dihedral_helix.h
  rm ../dihedral_hybrid.h
  rm ../dihedral_multi_harmonic.h
  rm ../dihedral_opls.h
  rm ../dump_bond.h
  rm ../fix_bond_break.h
  rm ../fix_bond_create.h
  rm ../fix_bond_swap.h
#  rm ../improper.h
  rm ../improper_cvff.h
  rm ../improper_harmonic.h
  rm ../improper_hybrid.h
  rm ../pair_lj_charmm_coul_charmm.h
  rm ../pair_lj_charmm_coul_charmm_implicit.h

endif
