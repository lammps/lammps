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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_nve_asphere_noforce.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "group.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEAsphereNoforce::FixNVEAsphereNoforce(LAMMPS *lmp, int narg, char **arg) :
  FixNVENoforce(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix nve/asphere/noforce command");

  time_integrate = 1;

  // error check

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!atom->ellipsoid_flag) 
    error->all(FLERR,"Fix nve/asphere/noforce requires atom style ellipsoid");
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereNoforce::init()
{
  FixNVENoforce::init();
  dtq = 0.5 * dtv;

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
	error->one(FLERR,"Fix nve/asphere/noforce requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereNoforce::initial_integrate(int vflag)
{
  AtomVecEllipsoid::Bonus *bonus;
  if (avec) bonus = avec->bonus;
  double **x = atom->x;
  double **v = atom->v;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *ellipsoid = atom->ellipsoid;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double *shape,*quat;
  double inertia[3],omega[3];

  // update positions and quaternions for all particles

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // principal moments of inertia

      shape = bonus[ellipsoid[i]].shape;
      quat = bonus[ellipsoid[i]].quat;
	
      inertia[0] = rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]) / 5.0;
      inertia[1] = rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]) / 5.0;
      inertia[2] = rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]) / 5.0;
      
      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion
      
      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);
    }
  }
}
