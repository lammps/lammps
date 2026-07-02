// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   AIREBO-BC: bond-centric modification of the P_ij term of AIREBO.
   Reference: J. Hur and S. J. Stuart, "Modified reactive empirical
   bond-order potential for heterogeneous bonding environments",
   J. Chem. Phys. 137, 054102 (2012).  https://doi.org/10.1063/1.4738879

   The bond-centric P_CC machinery lives in the PairAIREBO base class and is
   activated by bcflag; this variant only sets that flag.
------------------------------------------------------------------------- */

#include "pair_airebo_bc.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairAIREBObc::PairAIREBObc(LAMMPS *lmp) : PairAIREBO(lmp)
{
  bcflag = 1;
}
