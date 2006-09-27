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

#include "string.h"
#include "fix_energy.h"
#include "atom.h"
#include "neighbor.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

FixEnergy::FixEnergy(int narg, char **arg) : Fix(narg, arg)
{
  if (narg != 3) error->all("Illegal fix energy command");

  neigh_half_once = 1;
  nmax = 0;
  energy = NULL;
}

/* ---------------------------------------------------------------------- */

FixEnergy::~FixEnergy()
{
  memory->sfree(energy);
}

/* ---------------------------------------------------------------------- */

int FixEnergy::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnergy::init()
{
  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all("Pair style does not support dumping per-atom energy");

  eamstyle = 0;
  if (force->pair_match("eam")) eamstyle = 1;
  else if (force->pair_match("eam/alloy")) eamstyle = 1;
  else if (force->pair_match("eam/fs")) eamstyle = 1;

  // set communication size in comm class

  comm->maxreverse_fix = MAX(comm->maxreverse_fix,1);

  // warn if more than one energy fix
  
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"ENERGY") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one dump custom with an energy attribute");
}

/* ---------------------------------------------------------------------- */

void FixEnergy::dump()
{
  int i,j,k,n,itype,jtype,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double factor_coul,factor_lj,e;
  int *neighs;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(energy);
    nmax = atom->nmax;
    energy = (double *) memory->smalloc(nmax*sizeof(double),"energy:energy");
  }

  // clear energy array
  // n includes ghosts only if newton_pair flag is set

  if (force->newton_pair) n = atom->nlocal + atom->nghost;
  else n = atom->nlocal;

  for (i = 0; i < n; i++) energy[i] = 0.0;

  // if needed, build a half neighbor list

  if (!neighbor->half_every) neighbor->build_half();

  // compute pairwise energy for all atoms via pair->single()
  // use half neighbor list

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double **cutsq = force->pair->cutsq;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  Pair::One one;

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];
    
    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

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
	e = one.eng_coul + one.eng_vdwl;
	energy[i] += e;
	energy[j] += e;
      }
    }
  }

  // communicate energy between neighbor procs

  if (force->newton_pair) comm->reverse_comm_fix(this);

  // remove double counting of per-atom energy

  for (i = 0; i < nlocal; i++) energy[i] *= 0.5;

  // for EAM, include embedding function contribution to energy

  if (eamstyle) {
    int *type = atom->type;
    double fptmp,etmp;

    for (i = 0; i < nlocal; i++) {
      force->pair->single_embed(i,type[i],fptmp,1,etmp);
      energy[i] += etmp;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixEnergy::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixEnergy::unpack_reverse_comm(int n, int *list, double *buf)
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

int FixEnergy::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
