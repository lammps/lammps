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

#include <cstring>
#include "compute_meso_t_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoTAtom::ComputeMesoTAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Number of arguments for compute meso/t/atom command != 3");
  if ((atom->e_flag != 1) || (atom->cv_flag != 1))
          error->all(FLERR,"Compute meso/t/atom command requires atom_style with both energy and heat capacity (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  tvector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoTAtom::~ComputeMesoTAtom()
{
  memory->sfree(tvector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoTAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"meso/t/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute meso/t/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoTAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow tvector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(tvector);
    nmax = atom->nmax;
    tvector = (double *) memory->smalloc(nmax*sizeof(double),"tvector/atom:tvector");
    vector_atom = tvector;
  }

  double *e = atom->e;
  double *cv = atom->cv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              if (cv[i] > 0.0) {
                      tvector[i] = e[i] / cv[i];
              }
      }
      else {
              tvector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoTAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
