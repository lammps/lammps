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

#include "compute_sph_t_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPHTAtom::ComputeSPHTAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR,"Number of arguments for compute sph/t/atom command != 3");
  if ((atom->esph_flag != 1) || (atom->cv_flag != 1))
    error->all(FLERR,"Compute sph/t/atom command requires atom_style sph");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  tvector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSPHTAtom::~ComputeSPHTAtom()
{
  memory->sfree(tvector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPHTAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"meso/t/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute meso/t/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeSPHTAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow tvector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(tvector);
    nmax = atom->nmax;
    tvector = (double *) memory->smalloc(nmax*sizeof(double),"tvector/atom:tvector");
    vector_atom = tvector;
  }

  double *esph = atom->esph;
  double *cv = atom->cv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              if (cv[i] > 0.0) {
                      tvector[i] = esph[i] / cv[i];
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

double ComputeSPHTAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
