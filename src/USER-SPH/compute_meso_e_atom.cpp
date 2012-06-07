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
#include "compute_meso_e_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoEAtom::ComputeMesoEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Number of arguments for compute meso_e/atom command != 3");
  if (atom->e_flag != 1) error->all(FLERR,"compute meso_e/atom command requires atom_style with energy (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  evector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoEAtom::~ComputeMesoEAtom()
{
  memory->sfree(evector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoEAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"evector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute evector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoEAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow evector array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(evector);
    nmax = atom->nmax;
    evector = (double *) memory->smalloc(nmax*sizeof(double),"evector/atom:evector");
    vector_atom = evector;
  }

  double *e = atom->e;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              evector[i] = e[i];
      }
      else {
              evector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
