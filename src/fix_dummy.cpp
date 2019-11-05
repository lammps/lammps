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

#include "fix_dummy.h"
#include <cstring>
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDummy::FixDummy(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // process optional args
  // customize here and in setmask() by adding a new keyword from fix.h
  // only necessary if both of these are true:
  // (a) the real fix you are placeholding for defines the method
  // (b) the real fix will be defined so late in run initialization
  //     that the dummy fix will have already been processed by Modify::init()
  //     to add its index to its lists of fixes to invoke during timestepping

  initial_integrate_flag = final_integrate_flag = 0;
  pre_exchange_flag = pre_neighbor_flag = 0;
  pre_force_flag = post_force_flag = 0;
  end_of_step_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"initial_integrate") == 0) initial_integrate_flag = 1;
    else if (strcmp(arg[iarg],"final_integrate") == 0) final_integrate_flag = 1;
    else if (strcmp(arg[iarg],"final_integrate") == 0) final_integrate_flag = 1;
    else if (strcmp(arg[iarg],"final_integrate") == 0) final_integrate_flag = 1;
    else if (strcmp(arg[iarg],"final_integrate") == 0) final_integrate_flag = 1;
    else if (strcmp(arg[iarg],"final_integrate") == 0) final_integrate_flag = 1;
    else error->all(FLERR,"Illegal fix DUMMY command");
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

int FixDummy::setmask()
{
  int mask = 0;
  if (initial_integrate_flag) mask |= INITIAL_INTEGRATE;
  if (final_integrate_flag) mask |= FINAL_INTEGRATE;
  if (pre_exchange_flag) mask |= PRE_EXCHANGE;
  if (pre_neighbor_flag) mask |= PRE_NEIGHBOR;
  if (pre_force_flag) mask |= PRE_FORCE;
  if (post_force_flag) mask |= POST_FORCE;
  if (end_of_step_flag) mask |= END_OF_STEP;
  return mask;
}
