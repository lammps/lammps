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
#include "stdlib.h"
#include "string.h"
#include "fix_drag.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDrag::FixDrag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix drag command");

  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  xflag = yflag = zflag = 1;

  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else xc = force->numeric(FLERR,arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else yc = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else zc = force->numeric(FLERR,arg[5]);

  f_mag = force->numeric(FLERR,arg[6]);
  delta = force->numeric(FLERR,arg[7]);

  force_flag = 0;
  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixDrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDrag::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixDrag::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixDrag::post_force(int vflag)
{
  // apply drag force to atoms in group of magnitude f_mag
  // apply in direction (r-r0) if atom is further than delta away

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
  force_flag = 0;

  double dx,dy,dz,r,prefactor,fx,fy,fz;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dx = x[i][0] - xc;
      dy = x[i][1] - yc;
      dz = x[i][2] - zc;
      if (!xflag) dx = 0.0;
      if (!yflag) dy = 0.0;
      if (!zflag) dz = 0.0;
      domain->minimum_image(dx,dy,dz);
      r = sqrt(dx*dx + dy*dy + dz*dz);
      if (r > delta) {
        prefactor = f_mag/r;
        fx = prefactor*dx;
        fy = prefactor*dy;
        fz = prefactor*dz;
        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;
        ftotal[0] -= fx;
        ftotal[1] -= fy;
        ftotal[2] -= fz;
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixDrag::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total drag force on fix group
------------------------------------------------------------------------- */

double FixDrag::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(ftotal,ftotal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return ftotal_all[n];
}
