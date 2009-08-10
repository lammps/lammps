/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef AngleInclude
#include "angle_charmm.h"
#include "angle_cosine.h"
#include "angle_cosine_delta.h"
#include "angle_cosine_squared.h"
#include "angle_harmonic.h"
#include "angle_hybrid.h"
#include "angle_table.h"
#endif

#ifdef AngleClass
AngleStyle(charmm,AngleCharmm)
AngleStyle(cosine,AngleCosine)
AngleStyle(cosine/delta,AngleCosineDelta)
AngleStyle(cosine/squared,AngleCosineSquared)
AngleStyle(harmonic,AngleHarmonic)
AngleStyle(hybrid,AngleHybrid)
AngleStyle(table,AngleTable)
#endif

#ifdef AtomInclude
#include "atom_vec_angle.h"
#include "atom_vec_bond.h"
#include "atom_vec_full.h"
#include "atom_vec_molecular.h"
#endif

#ifdef AtomClass
AtomStyle(angle,AtomVecAngle)
AtomStyle(bond,AtomVecBond)
AtomStyle(full,AtomVecFull)
AtomStyle(molecular,AtomVecMolecular)
#endif

#ifdef BondInclude
#include "bond_fene.h"
#include "bond_fene_expand.h"
#include "bond_harmonic.h"
#include "bond_hybrid.h"
#include "bond_morse.h"
#include "bond_nonlinear.h"
#include "bond_quartic.h"
#include "bond_table.h"
#endif

#ifdef BondClass
BondStyle(fene,BondFENE)
BondStyle(fene/expand,BondFENEExpand)
BondStyle(harmonic,BondHarmonic)
BondStyle(hybrid,BondHybrid)
BondStyle(morse,BondMorse)
BondStyle(nonlinear,BondNonlinear)
BondStyle(quartic,BondQuartic)
BondStyle(table,BondTable)
#endif

#ifdef DihedralInclude
#include "dihedral_charmm.h"
#include "dihedral_harmonic.h"
#include "dihedral_helix.h"
#include "dihedral_hybrid.h"
#include "dihedral_multi_harmonic.h"
#include "dihedral_opls.h"
#endif

#ifdef DihedralClass
DihedralStyle(charmm,DihedralCharmm)
DihedralStyle(harmonic,DihedralHarmonic)
DihedralStyle(helix,DihedralHelix)
DihedralStyle(hybrid,DihedralHybrid)
DihedralStyle(multi/harmonic,DihedralMultiHarmonic)
DihedralStyle(opls,DihedralOPLS)
#endif

#ifdef DumpInclude
#include "dump_bond.h"
#endif

#ifdef DumpClass
DumpStyle(bond,DumpBond)
#endif

#ifdef FixInclude
#include "fix_bond_break.h"
#include "fix_bond_create.h"
#include "fix_bond_swap.h"
#endif

#ifdef FixClass
FixStyle(bond/break,FixBondBreak)
FixStyle(bond/create,FixBondCreate)
FixStyle(bond/swap,FixBondSwap)
#endif

#ifdef ImproperInclude
#include "improper_cvff.h"
#include "improper_harmonic.h"
#include "improper_hybrid.h"
#endif

#ifdef ImproperClass
ImproperStyle(cvff,ImproperCvff)
ImproperStyle(harmonic,ImproperHarmonic)
ImproperStyle(hybrid,ImproperHybrid)
#endif

#ifdef PairInclude
#include "pair_lj_charmm_coul_charmm.h"
#include "pair_lj_charmm_coul_charmm_implicit.h"
#endif

#ifdef PairClass
PairStyle(lj/charmm/coul/charmm,PairLJCharmmCoulCharmm)
PairStyle(lj/charmm/coul/charmm/implicit,PairLJCharmmCoulCharmmImplicit)
#endif
