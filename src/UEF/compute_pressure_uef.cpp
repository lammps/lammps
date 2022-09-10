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

#include "compute_pressure_uef.h"
#include <cstring>
#include "fix_nh_uef.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
 * Default values for the ext flags
 * ----------------------------------------------------------------------*/
ComputePressureUef::ComputePressureUef(LAMMPS *lmp, int narg, char **arg) :
  ComputePressure(lmp, narg, arg)
{
  ext_flags[0] = true;
  ext_flags[1] = true;
  ext_flags[2] = true;
  in_fix=false;
}

/* ----------------------------------------------------------------------
 *  Check for the uef fix
 * ----------------------------------------------------------------------*/
void ComputePressureUef::init()
{
  ComputePressure::init();
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
    error->all(FLERR,"Can't use compute pressure/uef without defining a fix nvt/npt/uef");
  ifix_uef=i;
  (dynamic_cast<FixNHUef*>(modify->fix[ifix_uef]))->get_ext_flags(ext_flags);

  if (strcmp(temperature->style,"temp/uef") != 0)
    error->warning(FLERR,"The temperature used in compute pressure/ued is not of style temp/uef");
}

/* ----------------------------------------------------------------------
 *  Compute pressure in the directions i corresponding to ext_flag[i]=true
 * ----------------------------------------------------------------------*/
double ComputePressureUef::compute_scalar()
{

  temperature->compute_scalar();
// if all pressures are external the scalar is found as normal
  if (ext_flags[0] && ext_flags[1] && ext_flags[2])
    return ComputePressure::compute_scalar();

// otherwise compute the full tensor and average desired components
  compute_vector();
  addstep(update->ntimestep+1);

  int k =0;
  scalar = 0.0;
  if (ext_flags[0]) {
    scalar += vector[0];
    k++;
  }
  if (ext_flags[1]) {
    scalar += vector[1];
    k++;
  }
  if (ext_flags[2]) {
    scalar += vector[2];
    k++;
  }

  if (k > 1) scalar /= k;
  return scalar;
}

/* ----------------------------------------------------------------------
   Compute the pressure tensor in the rotated coordinate system
------------------------------------------------------------------------- */
void ComputePressureUef::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' for "
               "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double *ke_tensor;
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    ke_tensor = temperature->vector;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (in_fix)
      virial_rot(virial,rot);
    else
    {
      double r[3][3];
      ( dynamic_cast<FixNHUef*>(modify->fix[ifix_uef]))->get_rot(r);
      virial_rot(virial,r);
    }
    if (keflag) {
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
    } else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
 * get the current rotation matrix and store it
------------------------------------------------------------------------- */
void ComputePressureUef::update_rot()
{
    ( dynamic_cast<FixNHUef*>(modify->fix[ifix_uef]))->get_rot(rot);
}

/* ----------------------------------------------------------------------
   Transform the pressure tensor to the rotated coordinate system
   [P]rot = Q.[P].Q^t
------------------------------------------------------------------------- */
void ComputePressureUef::virial_rot(double *x, const double r[3][3])
{

  double t[3][3];

  // [00 10 20 ] [ 0 3 4 ] [00 01 02 ]
  // [01 11 21 ] [ 3 1 5 ] [10 11 12 ]
  // [02 12 22 ] [ 4 5 2 ] [20 21 22 ]

  for (int k = 0; k<3; ++k)
  {
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
