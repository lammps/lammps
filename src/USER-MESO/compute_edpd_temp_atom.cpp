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
#include "compute_edpd_temp_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEDPDTempAtom::ComputeEDPDTempAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Number of arguments for compute edpd/temp/atom command != 3");
  if (atom->edpd_flag != 1) error->all(FLERR,"compute edpd/temp/atom command requires atom_style with temperature (e.g. edpd)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  temp_vector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeEDPDTempAtom::~ComputeEDPDTempAtom()
{
  memory->sfree(temp_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeEDPDTempAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"temp_vector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute temp_vector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeEDPDTempAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow temp_vector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(temp_vector);
    nmax = atom->nmax;
    temp_vector = (double *) memory->smalloc(nmax*sizeof(double),"temp_vector/atom:temp_vector");
    vector_atom = temp_vector;
  }

  double *edpd_temp = atom->edpd_temp;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              temp_vector[i] = edpd_temp[i];
      }
      else {
              temp_vector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEDPDTempAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
