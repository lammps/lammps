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
   Contributing author: Craig Tenney (UND) added support
                        for swapping atoms of different masses
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_viscosity.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// needs to be big, but not so big that lose precision when subtract velocity

#define BIG 1.0e10

/* ---------------------------------------------------------------------- */

FixViscosity::FixViscosity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix viscosity command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix viscosity command");

  scalar_flag = 1;
  global_freq = nevery;
  extscalar = 0;

  if (strcmp(arg[4],"x") == 0) vdim = 0;
  else if (strcmp(arg[4],"y") == 0) vdim = 1;
  else if (strcmp(arg[4],"z") == 0) vdim = 2;
  else error->all(FLERR,"Illegal fix viscosity command");

  if (strcmp(arg[5],"x") == 0) pdim = 0;
  else if (strcmp(arg[5],"y") == 0) pdim = 1;
  else if (strcmp(arg[5],"z") == 0) pdim = 2;
  else error->all(FLERR,"Illegal fix viscosity command");

  nbin = force->inumeric(FLERR,arg[6]);
  if (nbin % 2 || nbin <= 2) error->all(FLERR,"Illegal fix viscosity command");

  // optional keywords

  nswap = 1;
  vtarget = BIG;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"swap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix viscosity command");
      nswap = force->inumeric(FLERR,arg[iarg+1]);
      if (nswap <= 0) error->all(FLERR,"Fix viscosity swap value must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vtarget") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix viscosity command");
      if (strcmp(arg[iarg+1],"INF") == 0) vtarget = BIG;
      else vtarget = force->numeric(FLERR,arg[iarg+1]);
      if (vtarget <= 0.0)
        error->all(FLERR,"Fix viscosity vtarget value must be positive");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix viscosity command");
  }

  // initialize array sizes to nswap+1 so have space to shift values down

  pos_index = new int[nswap+1];
  neg_index = new int[nswap+1];
  pos_delta = new double[nswap+1];
  neg_delta = new double[nswap+1];

  p_exchange = 0.0;
}

/* ---------------------------------------------------------------------- */

FixViscosity::~FixViscosity()
{
  delete [] pos_index;
  delete [] neg_index;
  delete [] pos_delta;
  delete [] neg_delta;
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
  // warn if any fix ave/spatial comes after this fix
  // can cause glitch in averaging since ave will happen after swap

  int foundme = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i] == this) foundme = 1;
    if (foundme && strcmp(modify->fix[i]->style,"ave/spatial") == 0 && me == 0)
      error->warning(FLERR,"Fix viscosity comes before fix ave/spatial");
  }

  // set bounds of 2 slabs in pdim
  // only necessary for static box, else re-computed in end_of_step()
  // lo bin is always bottom bin
  // hi bin is just above half height

  if (domain->box_change == 0) {
    prd = domain->prd[pdim];
    boxlo = domain->boxlo[pdim];
    boxhi = domain->boxhi[pdim];
    double binsize = (boxhi-boxlo) / nbin;
    slablo_lo = boxlo;
    slablo_hi = boxlo + binsize;
    slabhi_lo = boxlo + (nbin/2)*binsize;
    slabhi_hi = boxlo + (nbin/2+1)*binsize;
  }

  periodicity = domain->periodicity[pdim];
}

/* ---------------------------------------------------------------------- */

void FixViscosity::end_of_step()
{
  int i,m,insert;
  double coord,delta;
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
    slabhi_lo = boxlo + (nbin/2)*binsize;
    slabhi_hi = boxlo + (nbin/2+1)*binsize;
  }

  // make 2 lists of up to nswap atoms with velocity closest to +/- vtarget
  // lists are sorted by closeness to vtarget
  // only consider atoms in the bottom/middle slabs
  // map atoms back into periodic box if necessary
  // insert = location in list to insert new atom

  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  npositive = nnegative = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      coord = x[i][pdim];
      if (coord < boxlo && periodicity) coord += prd;
      else if (coord >= boxhi && periodicity) coord -= prd;

      if (coord >= slablo_lo && coord < slablo_hi) {
        if (v[i][vdim] < 0.0) continue;
        delta = fabs(v[i][vdim] - vtarget);
        if (npositive < nswap || delta < pos_delta[nswap-1]) {
          for (insert = npositive-1; insert >= 0; insert--)
            if (delta > pos_delta[insert]) break;
          insert++;
          for (m = npositive-1; m >= insert; m--) {
            pos_delta[m+1] = pos_delta[m];
            pos_index[m+1] = pos_index[m];
          }
          pos_delta[insert] = delta;
          pos_index[insert] = i;
          if (npositive < nswap) npositive++;
        }
      }

      if (coord >= slabhi_lo && coord < slabhi_hi) {
        if (v[i][vdim] > 0.0) continue;
        delta = fabs(v[i][vdim] + vtarget);
        if (nnegative < nswap || delta < neg_delta[nswap-1]) {
          for (insert = nnegative-1; insert >= 0; insert--)
            if (delta > neg_delta[insert]) break;
          insert++;
          for (m = nnegative-1; m >= insert; m--) {
            neg_delta[m+1] = neg_delta[m];
            neg_index[m+1] = neg_index[m];
          }
          neg_delta[insert] = delta;
          neg_index[insert] = i;
          if (nnegative < nswap) nnegative++;
        }
      }
    }

  // loop over nswap pairs
  // find 2 global atoms with smallest delta in bottom/top slabs
  // BIG values are for procs with no atom to contribute
  // MINLOC also communicates which procs own them
  // exchange momenta between the 2 particles
  // if I own both particles just swap, else point2point comm of vel,mass

  double *mass = atom->mass;
  double *rmass = atom->rmass;

  int ipos,ineg;
  double sbuf[2],rbuf[2],vcm;

  double pswap = 0.0;
  mine[0].proc = mine[1].proc = me;
  int ipositive = 0;
  int inegative = 0;

  for (m = 0; m < nswap; m++) {
    if (ipositive < npositive) mine[0].value = pos_delta[ipositive];
    else mine[0].value = BIG;
    if (inegative < nnegative) mine[1].value = neg_delta[inegative];
    else mine[1].value = BIG;

    MPI_Allreduce(mine,all,2,MPI_DOUBLE_INT,MPI_MINLOC,world);

    if (all[0].value == BIG || all[1].value == BIG) continue;

    if (me == all[0].proc && me == all[1].proc) {
      ipos = pos_index[ipositive++];
      ineg = neg_index[inegative++];
      rbuf[0] = v[ipos][vdim];
      if (rmass) rbuf[1] = rmass[ipos];
      else rbuf[1] = mass[type[ipos]];
      sbuf[0] = v[ineg][vdim];
      if (rmass) sbuf[1] = rmass[ineg];
      else sbuf[1] = mass[type[ineg]];
      vcm = (sbuf[1]*sbuf[0] + rbuf[1]*rbuf[0]) / (sbuf[1] + rbuf[1]);
      v[ineg][vdim] = 2.0 * vcm - sbuf[0];
      v[ipos][vdim] = 2.0 * vcm - rbuf[0];
      pswap += rbuf[1] * (vcm - rbuf[0]) - sbuf[1] * (vcm - sbuf[0]);

    } else if (me == all[0].proc) {
      ipos = pos_index[ipositive++];
      sbuf[0] = v[ipos][vdim];
      if (rmass) sbuf[1] = rmass[ipos];
      else sbuf[1] = mass[type[ipos]];
      MPI_Sendrecv(sbuf,2,MPI_DOUBLE,all[1].proc,0,
                   rbuf,2,MPI_DOUBLE,all[1].proc,0,world,&status);
      vcm = (sbuf[1]*sbuf[0] + rbuf[1]*rbuf[0]) / (sbuf[1] + rbuf[1]);
      v[ipos][vdim] = 2.0 * vcm - sbuf[0];
      pswap += sbuf[1] * (vcm - sbuf[0]);

    } else if (me == all[1].proc) {
      ineg = neg_index[inegative++];
      sbuf[0] = v[ineg][vdim];
      if (rmass) sbuf[1] = rmass[ineg];
      else sbuf[1] = mass[type[ineg]];
      MPI_Sendrecv(sbuf,2,MPI_DOUBLE,all[0].proc,0,
                   rbuf,2,MPI_DOUBLE,all[0].proc,0,world,&status);
      vcm = (sbuf[1]*sbuf[0] + rbuf[1]*rbuf[0]) / (sbuf[1] + rbuf[1]);
      v[ineg][vdim] = 2.0 * vcm - sbuf[0];
      pswap -= sbuf[1] * (vcm - sbuf[0]);
    }
  }

  // tally momentum exchange from all swaps

  double pswap_all;
  MPI_Allreduce(&pswap,&pswap_all,1,MPI_DOUBLE,MPI_SUM,world);
  p_exchange += pswap_all;
}

/* ---------------------------------------------------------------------- */

double FixViscosity::compute_scalar()
{
  return p_exchange;
}
