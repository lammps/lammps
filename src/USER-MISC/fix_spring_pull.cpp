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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
   based on fix spring by: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_spring_pull.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

/* ---------------------------------------------------------------------- */

FixSpringPull::FixSpringPull(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix spring/pull command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 8;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  k_spring = force->numeric(FLERR,arg[3]);
  xflag = yflag = zflag = 1;
  if (strcmp(arg[4],"NULL") == 0) xflag = 0;
  else xc = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"NULL") == 0) yflag = 0;
  else yc = force->numeric(FLERR,arg[5]);
  if (strcmp(arg[6],"NULL") == 0) zflag = 0;
  else zc = force->numeric(FLERR,arg[6]);
  xv = force->numeric(FLERR,arg[7]);
  yv = force->numeric(FLERR,arg[8]);
  zv = force->numeric(FLERR,arg[9]);
  r0 = force->numeric(FLERR,arg[10]);
  if (r0 < 0) error->all(FLERR,"R0 < 0 for fix spring/pull command");

  ftotal[0] = ftotal[1] = ftotal[2] = ftotal[3] = 0.0;
  ftotal[4] = ftotal[5] = ftotal[6] = ftotal[7] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixSpringPull::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::init()
{
  masstotal = group->mass(igroup);
  
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::post_force(int vflag)
{
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  // fx,fy,fz = components of k * (r-r0) / masstotal
  
  double dx,dy,dz,fx,fy,fz,r,dr;

  // move reference point by one time step.
  xc += xv * update->dt;
  yc += yv * update->dt;
  zc += zv * update->dt;
  domain->minimum_image(xc,yc,zc);
  
  dx = xcm[0] - xc;
  dy = xcm[1] - yc;
  dz = xcm[2] - zc;
  domain->minimum_image(dx,dy,dz);
  
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  r = MAX(r,SMALL);
  dr = r - r0;

  fx = k_spring*dx*dr/r;
  fy = k_spring*dy*dr/r;
  fz = k_spring*dz*dr/r;
  ftotal[0] = -fx;
  ftotal[1] = -fy;
  ftotal[2] = -fz;
  ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
  ftotal[4] = r;
  ftotal[5] = xc;
  ftotal[6] = yc;
  ftotal[7] = zc;
  
  if (dr < 0.0) ftotal[3] = -ftotal[3]; 
  espring = 0.5*k_spring * dr*dr;

  fx /= masstotal;
  fy /= masstotal;
  fz /= masstotal;

  // apply restoring force to atoms in group

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  
  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringPull::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixSpringPull::compute_scalar()
{
  return espring;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixSpringPull::compute_vector(int n)
{
  return ftotal[n];
}
