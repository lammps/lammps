// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#include "compute_temp_uef.h"
#include <cstring>
#include "fix_nh_uef.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
 * Base constructor initialized to use rotation matrix
 * ----------------------------------------------------------------------*/
ComputeTempUef::ComputeTempUef(LAMMPS *lmp, int narg, char **arg) :
  ComputeTemp(lmp, narg, arg)
{
  rot_flag=true;
}

/* ----------------------------------------------------------------------
 *  Check for the uef fix
 * ----------------------------------------------------------------------*/
void ComputeTempUef::init()
{
  ComputeTemp::init();
  // check to make sure the other uef fix is on
  // borrowed from Pieter's nvt/sllod code
  int i=0;
  for (i=0; i<modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"nvt/uef")==0)
      break;
    if (strcmp(modify->fix[i]->style,"npt/uef")==0)
      break;
  }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute temp/uef without defining a fix nvt/npt/uef");
  ifix_uef=i;
}


/* ----------------------------------------------------------------------
   Compute the ke tensor in the proper coordinate system
------------------------------------------------------------------------- */
void ComputeTempUef::compute_vector()
{
  ComputeTemp::compute_vector();
  if (rot_flag) {
    double rot[3][3];
    ( (FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
    virial_rot(vector,rot);
  }

}

/* ----------------------------------------------------------------------
 * turn the rotation matrix on or off to properly account for the
 * coordinate system of the velocities
------------------------------------------------------------------------- */
void ComputeTempUef::yes_rot()
{
  rot_flag =true;
}
void ComputeTempUef::no_rot()
{
  rot_flag =false;
}

/* ----------------------------------------------------------------------
   Transform the pressure tensor to the rotated coordinate system
   [P]rot = Q.[P].Q^t
------------------------------------------------------------------------- */
void ComputeTempUef::virial_rot(double *x, const double r[3][3])
{

  double t[3][3];
  // [00 10 20 ] [ 0 3 4 ] [00 01 02 ]
  // [01 11 21 ] [ 3 1 5 ] [10 11 12 ]
  // [02 12 22 ] [ 4 5 2 ] [20 21 22 ]
  for (int k = 0; k<3; ++k) {
    t[0][k] = x[0]*r[0][k] + x[3]*r[1][k] + x[4]*r[2][k];
    t[1][k] = x[3]*r[0][k] + x[1]*r[1][k] + x[5]*r[2][k];
    t[2][k] = x[4]*r[0][k] + x[5]*r[1][k] + x[2]*r[2][k];
  }
  x[0] = r[0][0]*t[0][0] + r[1][0]*t[1][0] + r[2][0]*t[2][0];
  x[3] = r[0][0]*t[0][1] + r[1][0]*t[1][1] + r[2][0]*t[2][1];
  x[4] = r[0][0]*t[0][2] + r[1][0]*t[1][2] + r[2][0]*t[2][2];
  x[1] = r[0][1]*t[0][1] + r[1][1]*t[1][1] + r[2][1]*t[2][1];
  x[5] = r[0][1]*t[0][2] + r[1][1]*t[1][2] + r[2][1]*t[2][2];
  x[2] = r[0][2]*t[0][2] + r[1][2]*t[1][2] + r[2][2]*t[2][2];
}
