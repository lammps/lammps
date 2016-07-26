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

#include "np_granular_bin_newton_onesided.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighPairGranularBinNewtonOnesided::NeighPairGranularBinNewtonOnesided(LAMMPS *lmp) : 
  NeighPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NeighPairGranularBinNewtonOnesided::build(NeighList *list)
{
}
