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

#include "mpi.h"
#include "string.h"
#include "compute_group_group.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGroupGroup::ComputeGroupGroup(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute group/group command");

  scalar_flag = vector_flag = 1;
  size_vector = 3;
  extscalar = 1;
  extvector = 1;

  int n = strlen(arg[3]) + 1;
  group2 = new char[n];
  strcpy(group2,arg[3]);

  jgroup = group->find(group2);
  if (jgroup == -1) error->all("Compute group/group group ID does not exist");
  jgroupbit = group->bitmask[jgroup];

  vector = new double[3];
}

/* ---------------------------------------------------------------------- */

ComputeGroupGroup::~ComputeGroupGroup()
{
  delete [] group2;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::init()
{
  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all("Pair style does not support compute group/group");

  pair = force->pair;
  cutsq = force->pair->cutsq;

  // recheck that group 2 has not been deleted

  jgroup = group->find(group2);
  if (jgroup == -1) error->all("Compute group/group group ID does not exist");
  jgroupbit = group->bitmask[jgroup];

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeGroupGroup::compute_scalar()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  interact();
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::compute_vector()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  interact();
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::interact()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if i,j are not in 2 groups

  double one[4],all[4];
  one[0] = one[1] = one[2] = one[3] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) othergroupbit = jgroupbit;
    else if (mask[i] & jgroupbit) othergroupbit = groupbit;
    else continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      if (!(mask[j] & othergroupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

	// energy only computed once so tally full amount
	// force tally is jgroup acting on igroup

	if (newton_pair || j < nlocal) {
	  one[0] += eng;
	  if (othergroupbit == jgroupbit) {
	    one[1] += delx*fpair;
	    one[2] += dely*fpair;
	    one[3] += delz*fpair;
	  }
	  if (othergroupbit == groupbit) {
	    one[1] -= delx*fpair;
	    one[2] -= dely*fpair;
	    one[3] -= delz*fpair;
	  }

	// energy computed twice so tally half amount
	// only tally force if I own igroup atom

	} else {
	  one[0] += 0.5*eng;
	  if (othergroupbit == jgroupbit) {
	    one[1] += delx*fpair;
	    one[2] += dely*fpair;
	    one[3] += delz*fpair;
	  }
	}
      }
    }
  }

  MPI_Allreduce(one,all,4,MPI_DOUBLE,MPI_SUM,world);
  scalar = all[0];
  vector[0] = all[1]; vector[1] = all[2]; vector[2] = all[3];
}
