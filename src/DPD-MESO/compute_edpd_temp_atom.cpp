/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_edpd_temp_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEDPDTempAtom::ComputeEDPDTempAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  if (narg != 3)
    error->all(FLERR, 2, "Number of arguments for compute edpd/temp/atom command != 3");
  if (atom->edpd_flag != 1)
    error->all(FLERR, 2, "compute edpd/temp/atom command requires atom_style edpd");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  temp_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeEDPDTempAtom::~ComputeEDPDTempAtom()
{
  memory->sfree(temp_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeEDPDTempAtom::init()
{

  if ((comm->me == 0) && (modify->get_compute_by_style("^edpd/temp/atom").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeEDPDTempAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow temp_vector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(temp_vector);
    nmax = atom->nmax;
    temp_vector = (double *) memory->smalloc(nmax * sizeof(double), "edpd/temp/atom:temp_vector");
    vector_atom = temp_vector;
  }

  double *edpd_temp = atom->edpd_temp;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      temp_vector[i] = edpd_temp[i];
    } else {
      temp_vector[i] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEDPDTempAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
