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


#include "pointers.h"
#include "imbalance_group.h"

using namespace LAMMPS_NS;

int ImbalanceGroup::options(LAMMPS *lmp, int narg, char **arg)
{
  return 0;
}
 
void ImbalanceGroup::compute(LAMMPS *lmp, double *weight)
{
}
