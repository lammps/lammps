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

#ifdef AtomInclude
#include "atom_vec_peri.h"
#endif

#ifdef AtomClass
AtomStyle(peri,AtomVecPeri)
# endif

#ifdef FixInclude
#include "fix_peri_neigh.h"
#endif

#ifdef FixClass
FixStyle(PERI_NEIGH,FixPeriNeigh)
#endif

#ifdef PairInclude
#include "pair_peri_pmb.h"
#endif

#ifdef PairClass
PairStyle(peri/pmb,PairPeriPMB)
#endif

#ifdef ComputeInclude
#include "compute_damage_atom.h"
#endif

#ifdef ComputeClass
ComputeStyle(damage/atom,ComputeDamageAtom)
#endif

