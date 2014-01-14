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

#include "stdlib.h"
#include "string.h"
#include "fix_momentum.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

FixMomentum::FixMomentum(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix momentum command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix momentum command");

  linear = angular = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"linear") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix momentum command");
      linear = 1;
      xflag = force->inumeric(FLERR,arg[iarg+1]);
      yflag = force->inumeric(FLERR,arg[iarg+2]);
      zflag = force->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angular") == 0) {
      angular = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix momentum command");
  }

  if (linear == 0 && angular == 0)
    error->all(FLERR,"Illegal fix momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1) 
      error->all(FLERR,"Illegal fix momentum command");

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix momentum group has no atoms");
}

/* ---------------------------------------------------------------------- */

int FixMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMomentum::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixMomentum::end_of_step()
{
  if (linear) {
    double vcm[3];
    group->vcm(igroup,masstotal,vcm);

    // adjust velocities by vcm to zero linear momentum
    // only adjust a component if flag is set

    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xflag) v[i][0] -= vcm[0];
        if (yflag) v[i][1] -= vcm[1];
        if (zflag) v[i][2] -= vcm[2];
      }
  }

  if (angular) {
    double xcm[3],angmom[3],inertia[3][3],omega[3];
    group->xcm(igroup,masstotal,xcm);
    group->angmom(igroup,xcm,angmom);
    group->inertia(igroup,xcm,inertia);
    group->omega(angmom,inertia,omega);

    // adjust velocities to zero omega
    // vnew_i = v_i - w x r_i
    // must use unwrapped coords to compute r_i correctly

    double **x = atom->x;
    double **v = atom->v;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    double dx,dy,dz;
    double unwrap[3];

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - xcm[0];
        dy = unwrap[1] - xcm[1];
        dz = unwrap[2] - xcm[2];
        v[i][0] -= omega[1]*dz - omega[2]*dy;
        v[i][1] -= omega[2]*dx - omega[0]*dz;
        v[i][2] -= omega[0]*dy - omega[1]*dx;
      }
  }
}
