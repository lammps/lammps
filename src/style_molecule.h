/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef AngleInclude
#include "angle_charmm.h"
#include "angle_cosine.h"
#include "angle_cosine_squared.h"
#include "angle_harmonic.h"
#include "angle_hybrid.h"
#endif

#ifdef AngleClass
AngleStyle(charmm,AngleCharmm)
AngleStyle(cosine,AngleCosine)
AngleStyle(cosine/squared,AngleCosineSquared)
AngleStyle(harmonic,AngleHarmonic)
AngleStyle(hybrid,AngleHybrid)
#endif

#ifdef AtomInclude
#include "atom_angle.h"
#include "atom_bond.h"
#include "atom_full.h"
#include "atom_molecular.h"
#endif

#ifdef AtomClass
AtomStyle(angle,AtomAngle)
AtomStyle(bond,AtomBond)
AtomStyle(full,AtomFull)
AtomStyle(molecular,AtomMolecular)
#endif

#ifdef BondInclude
#include "bond_fene.h"
#include "bond_fene_expand.h"
#include "bond_harmonic.h"
#include "bond_hybrid.h"
#include "bond_morse.h"
#include "bond_nonlinear.h"
#include "bond_quartic.h"
#endif

#ifdef BondClass
BondStyle(fene,BondFENE)
BondStyle(fene/expand,BondFENEExpand)
BondStyle(harmonic,BondHarmonic)
BondStyle(hybrid,BondHybrid)
BondStyle(morse,BondMorse)
BondStyle(nonlinear,BondNonlinear)
BondStyle(quartic,BondQuartic)
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
#endif

#ifdef FixClass
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
