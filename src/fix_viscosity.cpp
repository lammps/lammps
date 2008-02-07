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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_viscosity.h"
#include "atom.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixViscosity::FixViscosity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all("Illegal fix viscosity command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix viscosity command");

  scalar_flag = 1;
  scalar_vector_freq = nevery;
  extscalar = 0;

  if (strcmp(arg[4],"x") == 0) vdim = 0;
  else if (strcmp(arg[4],"y") == 0) vdim = 1;
  else if (strcmp(arg[4],"z") == 0) vdim = 2;
  else error->all("Illegal fix viscosity command");

  if (strcmp(arg[5],"x") == 0) pdim = 0;
  else if (strcmp(arg[5],"y") == 0) pdim = 1;
  else if (strcmp(arg[5],"z") == 0) pdim = 2;
  else error->all("Illegal fix viscosity command");

  nbin = atoi(arg[6]);
  if (nbin < 3) error->all("Illegal fix viscosity command");

  flux = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixViscosity::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscosity::init()
{
  // set bounds of 2 slabs in pdim
  // only necessary for static box, else re-computed in end_of_step()
  // lo bin is always bottom bin
  // if nbin even, hi bin is just below half height
  // if nbin odd, hi bin straddles half height

  if (domain->box_change == 0) {
    prd = domain->prd[pdim];
    boxlo = domain->boxlo[pdim];
    boxhi = domain->boxhi[pdim];
    double binsize = (boxhi-boxlo) / nbin;
    slablo_lo = boxlo;
    slablo_hi = boxlo + binsize;
    slabhi_lo = boxlo + ((nbin-1)/2)*binsize;
    slabhi_hi = boxlo + ((nbin-1)/2 + 1)*binsize;
  }

  periodicity = domain->periodicity[pdim];
}

/* ---------------------------------------------------------------------- */

void FixViscosity::end_of_step()
{
  int i;
  double p,coord;
  MPI_Status status;
  struct {
    double value;
    int proc;
  } mine[2],all[2];

  // if box changes, recompute bounds of 2 slabs in pdim

  if (domain->box_change) {
    prd = domain->prd[pdim];
    boxlo = domain->boxlo[pdim];
    boxhi = domain->boxhi[pdim];
    double binsize = (boxhi-boxlo) / nbin;
    slablo_lo = boxlo;
    slablo_hi = boxlo + binsize;
    slabhi_lo = boxlo + ((nbin-1)/2)*binsize;
    slabhi_hi = boxlo + ((nbin-1)/2 + 1)*binsize;
  }

  // ipos,ineg = my 2 atoms with most pos/neg momenta in bottom/top slabs
  // map atom back into periodic box if necessary

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int ipos = -1;
  int ineg = -1;
  double pmin = BIG;
  double pmax = -BIG;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (mass) p = mass[type[i]] * v[i][vdim];
      else p = rmass[i] * v[i][vdim];
      coord = x[i][pdim];
      if (coord < boxlo && periodicity) coord += prd;
      else if (coord >= boxhi && periodicity) coord -= prd;
      if (coord >= slablo_lo && coord < slablo_hi) {
	if (p > pmax) {
	  pmax = p;
	  ipos = i;
	}
      }
      if (coord >= slabhi_lo && coord < slabhi_hi) {
	if (p < pmin) {
	  pmin = p;
	  ineg = i;
	}
      }
    }

  // find 2 global atoms with most pos/neg momenta in bottom/top slabs
  // MAXLOC also communicates which procs own them

  mine[0].value = pmax;
  mine[0].proc = me;
  mine[1].value = -pmin;
  mine[1].proc = me;

  MPI_Allreduce(mine,all,2,MPI_DOUBLE_INT,MPI_MAXLOC,world);

  // exchange momenta between the 2 particles
  // if I own both particles just swap, else point2point comm of mass,vel
  
  double sbuf[2],rbuf[2];

  if (me == all[0].proc && me == all[1].proc) {
    sbuf[0] = v[ineg][vdim];
    if (mass) sbuf[1] = mass[type[ineg]];
    else sbuf[1] = rmass[ineg];
    rbuf[0] = v[ipos][vdim];
    if (mass) rbuf[1] = mass[type[ipos]];
    else rbuf[1] = rmass[ipos];
    v[ineg][vdim] = rbuf[0] * rbuf[1]/sbuf[1];
    v[ipos][vdim] = sbuf[0] * sbuf[1]/rbuf[1];

  } else if (me == all[0].proc) {
    sbuf[0] = v[ipos][vdim];
    if (mass) sbuf[1] = mass[type[ipos]];
    else sbuf[1] = rmass[ipos];
    MPI_Sendrecv(sbuf,2,MPI_DOUBLE,all[1].proc,0,
		 rbuf,2,MPI_DOUBLE,all[1].proc,0,world,&status);
    v[ipos][vdim] = rbuf[0] * rbuf[1]/sbuf[1];

  } else if (me == all[1].proc) {
    sbuf[0] = v[ineg][vdim];
    if (mass) sbuf[1] = mass[type[ineg]];
    else sbuf[1] = rmass[ineg];
    MPI_Sendrecv(sbuf,2,MPI_DOUBLE,all[0].proc,0,
		 rbuf,2,MPI_DOUBLE,all[0].proc,0,world,&status);
    v[ineg][vdim] = rbuf[0] * rbuf[1]/sbuf[1];
  }

  // tally momentum flux
  // sign of all[1].value was flipped for MPI_Allreduce

  flux += all[0].value + all[1].value;
}

/* ---------------------------------------------------------------------- */

double FixViscosity::compute_scalar()
{
  return flux;
}
