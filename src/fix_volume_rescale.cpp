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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_volume_rescale.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "kspace.h"
#include "domain.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixVolRescale::FixVolRescale(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix volume/rescale command");

  box_change = 1;

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix volume/rescale command");

  xflag = yflag = zflag = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix volume/rescale command");
      if (domain->xperiodic == 0) 
	error->all("Cannot fix volume/rescale on a non-periodic boundary");
      xflag = 1;
      xlo_stop = atof(arg[iarg+1]);
      xhi_stop = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix volume/rescale command");
      if (domain->yperiodic == 0) 
	error->all("Cannot fix volume/rescale on a non-periodic boundary");
      yflag = 1;
      ylo_stop = atof(arg[iarg+1]);
      yhi_stop = atof(arg[iarg+2]);
      iarg += 3;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix volume/rescale command");
      if (domain->zperiodic == 0) 
	error->all("Cannot fix volume/rescale on a non-periodic boundary");
      zflag = 1;
      zlo_stop = atof(arg[iarg+1]);
      zhi_stop = atof(arg[iarg+2]);
      iarg += 3;
    } else error->all("Illegal fix volume/rescale command");
  }

  if (domain->triclinic)
    error->all("Cannot use fix volume/rescale with triclinic box");

  nrigid = 0;
  rfix = NULL;
}

/* ---------------------------------------------------------------------- */

FixVolRescale::~FixVolRescale()
{
  delete [] rfix;
}

/* ---------------------------------------------------------------------- */

int FixVolRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixVolRescale::init()
{
  xlo_start = domain->boxxlo;
  xhi_start = domain->boxxhi;
  ylo_start = domain->boxylo;
  yhi_start = domain->boxyhi;
  zlo_start = domain->boxzlo;
  zhi_start = domain->boxzhi;

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any fix rigid exist so rigid bodies can be rescaled
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	strcmp(modify->fix[i]->style,"poems") == 0) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	  strcmp(modify->fix[i]->style,"poems") == 0) rfix[nrigid++] = i;
  }
}

/* ---------------------------------------------------------------------- */

void FixVolRescale::end_of_step()
{
  int i;
  double oldlo,oldhi,newlo,newhi,ratio;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (xflag) {
    oldlo = domain->boxxlo;
    oldhi = domain->boxxhi;
    newlo = xlo_start + delta * (xlo_stop-xlo_start);
    newhi = xhi_start + delta * (xhi_stop-xhi_start);
    ratio = (newhi - newlo) / domain->xprd;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	x[i][0] = newlo + (x[i][0] - oldlo) * ratio;
    domain->boxxlo = newlo;
    domain->boxxhi = newhi;
    domain->xprd = newhi - newlo;
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->dilate(0,oldlo,oldhi,newlo,newhi);
  }

  if (yflag) {
    oldlo = domain->boxylo;
    oldhi = domain->boxyhi;
    newlo = ylo_start + delta * (ylo_stop-ylo_start);
    newhi = yhi_start + delta * (yhi_stop-yhi_start);
    ratio = (newhi - newlo) / domain->yprd;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	x[i][1] = newlo + (x[i][1] - oldlo) * ratio;
    domain->boxylo = newlo;
    domain->boxyhi = newhi;
    domain->yprd = newhi - newlo;
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->dilate(1,oldlo,oldhi,newlo,newhi);
  }

  if (zflag) {
    oldlo = domain->boxzlo;
    oldhi = domain->boxzhi;
    newlo = zlo_start + delta * (zlo_stop-zlo_start);
    newhi = zhi_start + delta * (zhi_stop-zhi_start);
    ratio = (newhi - newlo) / domain->zprd;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	x[i][2] = newlo + (x[i][2] - oldlo) * ratio;
    domain->boxzlo = newlo;
    domain->boxzhi = newhi;
    domain->zprd = newhi - newlo;
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	  modify->fix[rfix[i]]->dilate(2,oldlo,oldhi,newlo,newhi);
  }

  // redo KSpace coeffs since volume has changed

  if (kspace_flag) force->kspace->setup();
}
