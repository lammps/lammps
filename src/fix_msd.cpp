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
#include "fix_msd.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixMSD::FixMSD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix msd command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix msd command");
  first = 1;
  restart_peratom = 1;

  MPI_Comm_rank(world,&me);
  if (me == 0) {
    fp = fopen(arg[4],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix msd file %s",arg[4]);
      error->one(str);
    }
  }

  if (me == 0) {
    fprintf(fp,"# Mean-squared Displacement for group %s\n",
	    group->names[igroup]);
    fprintf(fp,"# TimeStep x y z total\n");
  }

  // optional args

  comflag = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix msd command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all("Illegal fix msd command");
      iarg += 2;
    } else error->all("Illegal fix msd command");
  }

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // cm = original center of mass

  double cm[3];

  if (comflag) {
    masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,cm);
  }

  // xoriginal = initial unwrapped positions of atoms
  // relative to center of mass if comflag is set

  double **x = atom->x;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],xoriginal[i]);
      if (comflag) {
	xoriginal[i][0] -= cm[0];
	xoriginal[i][1] -= cm[1];
	xoriginal[i][2] -= cm[2];
      }
    } else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }

  // nmsd = # of atoms in group

  nmsd = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) nmsd++;

  int nmsd_all;
  MPI_Allreduce(&nmsd,&nmsd_all,1,MPI_INT,MPI_SUM,world);
  nmsd = nmsd_all;

  if (nmsd == 0) error->all("Fix msd group has no atoms");
}

/* ---------------------------------------------------------------------- */

FixMSD::~FixMSD()
{
  if (me == 0) fclose(fp);

  // unregister callbacks to this fix from Atom class
 
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored array

  memory->destroy_2d_double_array(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixMSD::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMSD::init()
{
  // warn if more than one msd fix

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"msd") == 0) count++;
  if (count > 1 && me == 0) error->warning("More than one fix msd");
}

/* ---------------------------------------------------------------------- */

void FixMSD::setup(int vflag)
{
  if (first) end_of_step();
  first = 0;
}

/* ---------------------------------------------------------------------- */

void FixMSD::end_of_step()
{
  // cm = current center of mass

  double cm[3];

  if (comflag) {
    group->mass(igroup);
    group->xcm(igroup,masstotal,cm);
  }

  double **x = atom->x;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  int xbox,ybox,zbox;
  double dx,dy,dz;
  double msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;
  
  // dx,dy,dz = displacement of atom from original position
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	xbox = (image[i] & 1023) - 512;
	ybox = (image[i] >> 10 & 1023) - 512;
	zbox = (image[i] >> 20) - 512;
	if (comflag) {
	  dx = x[i][0] + xbox*xprd - cm[0] - xoriginal[i][0];
	  dy = x[i][1] + ybox*yprd - cm[1] - xoriginal[i][1];
	  dz = x[i][2] + zbox*zprd - cm[2] - xoriginal[i][2];
	} else {
	  dx = x[i][0] + xbox*xprd - xoriginal[i][0];
	  dy = x[i][1] + ybox*yprd - xoriginal[i][1];
	  dz = x[i][2] + zbox*zprd - xoriginal[i][2];
	}
	msd[0] += dx*dx;
	msd[1] += dy*dy;
	msd[2] += dz*dz;
	msd[3] += dx*dx + dy*dy + dz*dz;
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	xbox = (image[i] & 1023) - 512;
	ybox = (image[i] >> 10 & 1023) - 512;
	zbox = (image[i] >> 20) - 512;
	if (comflag) {
	  dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - 
	    cm[0] - xoriginal[i][0];
	  dy = x[i][1] + h[1]*ybox + h[3]*zbox - cm[1] - xoriginal[i][1];
	  dz = x[i][2] + h[2]*zbox - cm[2] - xoriginal[i][2];
	} else {
	  dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox - xoriginal[i][0];
	  dy = x[i][1] + h[1]*ybox + h[3]*zbox - xoriginal[i][1];
	  dz = x[i][2] + h[2]*zbox - xoriginal[i][2];
	}
	msd[0] += dx*dx;
	msd[1] += dy*dy;
	msd[2] += dz*dz;
	msd[3] += dx*dx + dy*dy + dz*dz;
      }
  }

  double msd_all[4];
  MPI_Allreduce(msd,msd_all,4,MPI_DOUBLE,MPI_SUM,world);
  msd_all[0] /= nmsd;
  msd_all[1] /= nmsd;
  msd_all[2] /= nmsd;
  msd_all[3] /= nmsd;

  if (me == 0) {
    fprintf(fp,"%d %g %g %g %g\n",update->ntimestep,
	    msd_all[0],msd_all[1],msd_all[2],msd_all[3]);
    fflush(fp);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMSD::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixMSD::grow_arrays(int nmax)
{
  xoriginal =
    memory->grow_2d_double_array(xoriginal,nmax,3,"msd:xoriginal");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixMSD::copy_arrays(int i, int j)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixMSD::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixMSD::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixMSD::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixMSD::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixMSD::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixMSD::size_restart(int nlocal)
{
  return 4;
}
