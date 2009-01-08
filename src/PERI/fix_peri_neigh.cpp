/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "fix_peri_neigh.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

#include "comm.h"
#include "update.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixPeriNeigh::FixPeriNeigh(LAMMPS *lmp,int narg, char **arg) : 
  Fix(lmp, narg, arg)
{
  restart_peratom = 1;
  first = 1;

  // perform initial allocation of atom-based arrays
  // register with atom class
  // set maxpartner = 1 as placeholder

  maxpartner = 1;
  npartner = NULL;
  partner = NULL;
  r0 = NULL;
  vinter = NULL;

  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so atom migration is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;
}

/* ---------------------------------------------------------------------- */

FixPeriNeigh::~FixPeriNeigh()
{
  // if atom class still exists:
  // unregister this fix so atom class doesn't invoke it any more

  if (atom) atom->delete_callback(id,0);
  if (atom) atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->sfree(npartner);
  memory->destroy_2d_int_array(partner);
  memory->destroy_2d_double_array(r0);
  memory->sfree(vinter);
}

/* ---------------------------------------------------------------------- */

int FixPeriNeigh::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPeriNeigh::init()
{
  // need a full neighbor list once

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0; 
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */
 
void FixPeriNeigh::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   create initial list of neighbor partners via call to neighbor->build()
   must be done in setup (not init) since fix init comes before neigh init
------------------------------------------------------------------------- */

void FixPeriNeigh::setup(int vflag)
{
  int i,j,ii,jj,itype,jtype,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh;
  int **firstneigh;

  double **x = atom->x;
  double *vfrac = atom->vfrac;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;

  // only build list of bonds on very first run

  if (!first) return;
  first = 0;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // scan neighbor list to set maxpartner

  Pair *anypair = force->pair_match("peri",0);
  double **cutsq = anypair->cutsq;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq <= cutsq[itype][jtype]) npartner[i]++;
    }
  }

  maxpartner = 0;
  for (i = 0; i < nlocal; i++) maxpartner = MAX(maxpartner,npartner[i]);
  int maxall;
  MPI_Allreduce(&maxpartner,&maxall,1,MPI_INT,MPI_MAX,world);
  maxpartner = maxall;

  // realloc arrays with correct value for maxpartner

  memory->destroy_2d_int_array(partner);
  memory->destroy_2d_double_array(r0);
  memory->sfree(npartner);

  npartner = NULL;
  partner = NULL;
  r0 = NULL;
  grow_arrays(atom->nmax);

  // create partner list and r0 values from neighbor list
  // compute vinter for each atom

  for (i = 0; i < nlocal; i++) {
    npartner[i] = 0;
    vinter[i] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
 
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq <= cutsq[itype][jtype]) {
	partner[i][npartner[i]] = tag[j];
	r0[i][npartner[i]] = sqrt(rsq);
	npartner[i]++;
        vinter[i] += vfrac[j];
      }
    }
  }

  // bond statistics

  int n = 0;
  for (i = 0; i < nlocal; i++) n += npartner[i];
  int nall;
  MPI_Allreduce(&n,&nall,1,MPI_INT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Peridynamic bonds:\n");
      fprintf(screen,"  total # of bonds = %d\n",nall);
      fprintf(screen,"  bonds/atom = %g\n",nall/atom->natoms);
    }
    if (logfile) {
      fprintf(logfile,"Peridynamic bonds:\n");
      fprintf(logfile,"  total # of bonds = %d\n",nall);
      fprintf(logfile,"  bonds/atom = %g\n",nall/atom->natoms);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPeriNeigh::memory_usage()
{
  int nmax = atom->nmax;
  int bytes = nmax * sizeof(int);
  bytes += nmax*maxpartner * sizeof(int);
  bytes += nmax*maxpartner * sizeof(double);
  bytes += nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPeriNeigh::grow_arrays(int nmax)
{
  npartner = (int *) memory->srealloc(npartner,nmax*sizeof(int),
				      "peri_neigh:npartner");
  partner = memory->grow_2d_int_array(partner,nmax,maxpartner,
				      "peri_neigh:partner");
  r0 = memory->grow_2d_double_array(r0,nmax,maxpartner,"peri_neigh:r0");
  vinter = (double *) memory->srealloc(vinter,nmax*sizeof(double),
				       "peri_neigh:vinter");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPeriNeigh::copy_arrays(int i, int j)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    r0[j][m] = r0[i][m];
  }
  vinter[j] = vinter[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPeriNeigh::pack_exchange(int i, double *buf)
{
  // compact list by eliminating partner = 0 entries
  // set buf[0] after compaction

  int m = 1;
  for (int n = 0; n < npartner[i]; n++) {
    if (partner[i][n] == 0) continue;
    buf[m++] = partner[i][n];
    buf[m++] = r0[i][n];
  }

  buf[0] = m/2;
  buf[m++] = vinter[i];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixPeriNeigh::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    r0[nlocal][n] = buf[m++];
  }
  vinter[nlocal] = buf[m++];
  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPeriNeigh::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = 2*npartner[i] + 2;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    buf[m++] = r0[i][n];
  }
  buf[m++] = vinter[i];
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPeriNeigh::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    r0[nlocal][n] = extra[nlocal][m++];
  }
  vinter[nlocal] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPeriNeigh::maxsize_restart()
{
  return 2*maxpartner + 3;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPeriNeigh::size_restart(int nlocal)
{
  return 2*npartner[nlocal] + 3;
}
