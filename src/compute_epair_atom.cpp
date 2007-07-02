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
#include "compute_epair_atom.h"
#include "atom.h"
#include "neighbor.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ComputeEpairAtom::ComputeEpairAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute epair/atom command");

  peratom_flag = 1;
  size_peratom = 0;
  comm_reverse = 1;
  neigh_half_once = 1;

  nmax = 0;
  energy = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEpairAtom::~ComputeEpairAtom()
{
  memory->sfree(energy);
}

/* ---------------------------------------------------------------------- */

void ComputeEpairAtom::init()
{
  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all("Pair style does not support computing per-atom energy");

  if (force->pair_match("eam")) eamstyle = 1;
  else eamstyle = 0;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"epair/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute epair/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEpairAtom::compute_peratom()
{
  int i,j,k,n,itype,jtype,numneigh,iflag;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double factor_coul,factor_lj,eng;
  int *neighs;
  Pair::One one;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(energy);
    nmax = atom->nmax;
    energy = (double *) 
      memory->smalloc(nmax*sizeof(double),"compute/epair/atom:energy");
    scalar_atom = energy;
  }

  // clear energy array
  // n includes ghosts only if newton_pair flag is set

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  if (newton_pair) n = nlocal + atom->nghost;
  else n = nlocal;

  for (i = 0; i < n; i++) energy[i] = 0.0;

  // if needed, build a half neighbor list

  if (!neighbor->half_every) neighbor->build_half();

  // compute pairwise energy for atoms via pair->single()
  // if neither atom is in compute group, skip that pair

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double **cutsq = force->pair->cutsq;

  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    iflag = mask[i] & groupbit;
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];
    
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];
      if (iflag == 0 && (mask[j] & groupbit) == 0) continue;

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	force->pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,1,one);
	eng = one.eng_coul + one.eng_vdwl;
	energy[i] += eng;
	if (newton_pair || j < nlocal) energy[j] += eng;
      }
    }
  }

  // communicate energy between neigchbor procs

  if (newton_pair) comm->reverse_comm_compute(this);

  // remove double counting of per-atom energy

  for (i = 0; i < nlocal; i++) energy[i] *= 0.5;

  // for EAM, include embedding function contribution to energy
  // only for atoms in compute group

  if (eamstyle) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	force->pair->single_embed(i,type[i],eng);
	energy[i] += eng;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeEpairAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeEpairAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeEpairAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
