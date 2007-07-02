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
#include "compute_ebond_atom.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "bond.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ComputeEbondAtom::ComputeEbondAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute ebond/atom command");

  peratom_flag = 1;
  size_peratom = 0;
  comm_reverse = 1;

  nmax = 0;
  energy = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEbondAtom::~ComputeEbondAtom()
{
  memory->sfree(energy);
}

/* ---------------------------------------------------------------------- */

void ComputeEbondAtom::init()
{
  if (force->bond == NULL)
    error->all("Bond style does not support computing per-atom energy");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"ebond/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute ebond/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEbondAtom::compute_peratom()
{
  int i,n,i1,i2,type;
  double delx,dely,delz,rsq,fforce,eng;

  // grow energy array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(energy);
    nmax = atom->nmax;
    energy = (double *) 
      memory->smalloc(nmax*sizeof(double),"compute/epair/atom:energy");
    scalar_atom = energy;
  }

  // clear energy array
  // n includes ghosts only if newton_bond flag is set

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  if (newton_bond) n = nlocal + atom->nghost;
  else n = nlocal;

  for (i = 0; i < n; i++) energy[i] = 0.0;

  // compute bond energy for atoms via bond->single()
  // if neither atom is in compute group, skip that bond

  double **x = atom->x;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    if ((mask[i1] & groupbit) == 0 && (mask[i2] & groupbit) == 0) continue;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);
      
    rsq = delx*delx + dely*dely + delz*delz;
    force->bond->single(type,rsq,i1,i2,0,fforce,eng);
    energy[i] += eng;
    if (newton_bond || i2 < nlocal) energy[i2] += eng;
  }

  // communicate energy between neigchbor procs

  if (newton_bond) comm->reverse_comm_compute(this);

  // remove double counting of per-atom energy

  for (i = 0; i < nlocal; i++) energy[i] *= 0.5;
}

/* ---------------------------------------------------------------------- */

int ComputeEbondAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeEbondAtom::unpack_reverse_comm(int n, int *list, double *buf)
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

int ComputeEbondAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
