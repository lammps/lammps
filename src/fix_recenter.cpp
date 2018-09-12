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
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_recenter.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "modify.h"
#include "comm.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{BOX,LATTICE,FRACTION};

/* ---------------------------------------------------------------------- */

FixRecenter::FixRecenter(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix recenter command");

  xcom = ycom = zcom = 0.0;
  xflag = yflag = zflag = 1;
  xinitflag = yinitflag = zinitflag = 0;
  shift[0] = shift[1] = shift[2] = 0.0;
  distance = 0.0;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  extscalar = 1;
  extvector = 1;
  global_freq = 1;
  dynamic_group_allow = 1;

  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else if (strcmp(arg[3],"INIT") == 0) xinitflag = 1;
  else xcom = force->numeric(FLERR,arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else if (strcmp(arg[4],"INIT") == 0) yinitflag = 1;
  else ycom = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else if (strcmp(arg[5],"INIT") == 0) zinitflag = 1;
  else zcom = force->numeric(FLERR,arg[5]);

  // optional args

  group2bit = groupbit;
  scaleflag = LATTICE;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"shift") == 0) {
      int igroup2 = group->find(arg[iarg+1]);
      if (igroup2 < 0) error->all(FLERR,"Could not find fix recenter group ID");
      group2bit = group->bitmask[igroup2];
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"fraction") == 0) scaleflag = FRACTION;
      else error->all(FLERR,"Illegal fix recenter command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix recenter command");
  }

  // scale xcom,ycom,zcom

  double xscale,yscale,zscale;
  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  xcom *= xscale;
  ycom *= yscale;
  zcom *= zscale;

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix recenter group has no atoms");
}

/* ---------------------------------------------------------------------- */

int FixRecenter::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRecenter::init()
{
  // warn if any integrate fix comes after this one

  int after = 0;
  int flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(id,modify->fix[i]->id) == 0) after = 1;
    else if ((modify->fmask[i] & INITIAL_INTEGRATE) && after) flag = 1;
  }
  if (flag && comm->me == 0)
    error->warning(FLERR,"Fix recenter should come after all other "
                   "integration fixes");

  masstotal = group->mass(igroup);

  // if any components of requested COM were INIT, store initial COM

  if (xinitflag || yinitflag || zinitflag) {
    double xcm[3];
    group->xcm(igroup,masstotal,xcm);
    xinit = xcm[0];
    yinit = xcm[1];
    zinit = xcm[2];
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRecenter::initial_integrate(int /*vflag*/)
{
  // target COM
  // bounding box around domain works for both orthogonal and triclinic

  double xtarget,ytarget,ztarget;
  double *bboxlo,*bboxhi;

  if (scaleflag == FRACTION) {
    if (domain->triclinic == 0) {
      bboxlo = domain->boxlo;
      bboxhi = domain->boxhi;
    } else {
      bboxlo = domain->boxlo_bound;
      bboxhi = domain->boxhi_bound;
    }
  }

  if (xinitflag) xtarget = xinit;
  else if (scaleflag == FRACTION)
    xtarget = bboxlo[0] + xcom*(bboxhi[0] - bboxlo[0]);
  else xtarget = xcom;

  if (yinitflag) ytarget = yinit;
  else if (scaleflag == FRACTION)
    ytarget = bboxlo[1] + ycom*(bboxhi[1] - bboxlo[1]);
  else ytarget = ycom;

  if (zinitflag) ztarget = zinit;
  else if (scaleflag == FRACTION)
    ztarget = bboxlo[2] + zcom*(bboxhi[2] - bboxlo[2]);
  else ztarget = zcom;

  // current COM

  double xcm[3];
  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);

  group->xcm(igroup,masstotal,xcm);

  // shift coords by difference between actual COM and requested COM

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  shift[0] = xflag ? (xtarget - xcm[0]) : 0.0;
  shift[1] = yflag ? (ytarget - xcm[1]) : 0.0;
  shift[2] = zflag ? (ztarget - xcm[2]) : 0.0;
  distance = sqrt(shift[0]*shift[0] + shift[1]*shift[1] + shift[2]*shift[2]);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & group2bit) {
      x[i][0] += shift[0];
      x[i][1] += shift[1];
      x[i][2] += shift[2];
    }
}

/* ---------------------------------------------------------------------- */

void FixRecenter::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  // outermost level - operate recenter
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

double FixRecenter::compute_scalar()
{
  return distance;
}

/* ---------------------------------------------------------------------- */

double FixRecenter::compute_vector(int n)
{
  return shift[n];
}
