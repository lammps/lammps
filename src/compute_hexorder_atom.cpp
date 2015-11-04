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
   Contributing author:  Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "compute_hexorder_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::ComputeHexOrderAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute hexorder/atom command");

  double cutoff = force->numeric(FLERR,arg[3]);
  cutsq = cutoff*cutoff;

  ncol = 2;

  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  q6array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeHexOrderAtom::~ComputeHexOrderAtom()
{
  memory->destroy(q6array);
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute hexorder/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute hexorder/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"hexorder/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute hexorder/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeHexOrderAtom::compute_peratom()
{
  int i,j,m,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(q6array);
    nmax = atom->nmax;
    memory->create(q6array,nmax,ncol,"hexorder/atom:q6array");
    array_atom = q6array;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter(s) for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    double* q6 = q6array[i];
    q6[0] = q6[1] = 0.0;
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      double usum = 0.0;
      double vsum = 0.0;
      int ncount = 0;
      
      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	j &= NEIGHMASK;
	
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	if (rsq < cutsq) {
	  double u, v;
	  calc_q6(delx, dely, u, v);
	  usum += u;
	  vsum += v;
	  ncount++;
	}
      }
      if (ncount > 0) {
	double ninv = 1.0/ncount ;
	q6[0] = usum*ninv;
	q6[1] = vsum*ninv;
      }
    }
  }
}

inline void ComputeHexOrderAtom::calc_q6(double delx, double dely, double &u, double &v) {
  double rinv = 1.0/sqrt(delx*delx+dely*dely);
  double x = delx*rinv;
  double y = dely*rinv;
  double a = x*x;
  double b1 = y*y;
  double b2 = b1*b1;
  double b3 = b2*b1;
  u = ((  a - 15*b1)*a + 15*b2)*a - b3;
  v = ((6*a - 20*b1)*a +  6*b2)*x*y;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeHexOrderAtom::memory_usage()
{
  double bytes = ncol*nmax * sizeof(double);
  return bytes;
}
