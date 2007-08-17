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
#include "compute_coord_atom.h"
#include "atom.h"
#include "modify.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::ComputeCoordAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute centro/atom command");

  cutoff = atof(arg[3]);

  peratom_flag = 1;
  size_peratom = 0;
  neigh_full_once = 1;

  nmax = 0;
  coordination = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::~ComputeCoordAtom()
{
  memory->sfree(coordination);
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::init()
{
  if (force->pair == NULL || cutoff > force->pair->cutforce) 
    error->all("Compute coord/atom cutoff is longer than pairwise cutoff");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"coord/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute coord/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::compute_peratom()
{
  int j,k,n,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighs;

  // grow coordination array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(coordination);
    nmax = atom->nmax;
    coordination = (double *) 
      memory->smalloc(nmax*sizeof(double),"compute/coord/atom:coordination");
    scalar_atom = coordination;
  }

  // if needed, build a full neighbor list

  if (!neighbor->full_every) neighbor->build_full();

  // compute coordination number for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double cutsq = cutoff*cutoff;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      neighs = neighbor->firstneigh_full[i];
      numneigh = neighbor->numneigh_full[i];

      n = 0;
      for (k = 0; k < numneigh; k++) {
	j = neighs[k];
	if (j >= nall) j %= nall;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	if (rsq < cutsq) n++;
      }

      coordination[i] = n;
    }
}
/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeCoordAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
