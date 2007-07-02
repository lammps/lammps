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
#include "compute_etotal_atom.h"
#include "atom.h"
#include "modify.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEtotalAtom::ComputeEtotalAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute etotal/atom command");

  // store epair ID used by energy computation

  int n = strlen(arg[3]) + 1;
  id_pre = new char[n];
  strcpy(id_pre,arg[3]);

  peratom_flag = 1;
  size_peratom = 0;

  nmax = 0;
  etotal = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEtotalAtom::~ComputeEtotalAtom()
{
  delete [] id_pre;
  memory->sfree(etotal);
}

/* ---------------------------------------------------------------------- */

void ComputeEtotalAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"etotal/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute etotal/atom");

  // set epair Compute used by this compute

  int icompute = modify->find_compute(id_pre);
  if (icompute < 0)
    error->all("Could not find compute etotal/atom pre-compute ID");
  compute_epair = modify->compute[icompute];

  if (groupbit != compute_epair->groupbit && comm->me == 0)
    error->warning("Group for compute etotal and its epair are not the same");
}

/* ---------------------------------------------------------------------- */

void ComputeEtotalAtom::compute_peratom()
{
  // grow etotal array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(etotal);
    nmax = atom->nmax;
    etotal = (double *) 
      memory->smalloc(nmax*sizeof(double),"compute/etotal/atom:etotal");
    scalar_atom = etotal;
  }

  // compute total energy for each atom in group

  double *epair = compute_epair->scalar_atom;
  double mvv2e = force->mvv2e;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double ke;

  if (mass)
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ke = 0.5 * mvv2e * mass[type[i]] *
	  (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
	etotal[i] = ke + epair[i];
      }
    }
  else
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	ke = 0.5 * mvv2e * rmass[i] *
	  (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
	etotal[i] = ke + epair[i];
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

int ComputeEtotalAtom::memory_usage()
{
  int bytes = nmax * sizeof(double);
  return bytes;
}
