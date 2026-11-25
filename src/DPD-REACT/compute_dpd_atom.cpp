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

/* ----------------------------------------------------------------------
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "compute_dpd_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDpdAtom::ComputeDpdAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), dpdAtom(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute dpd/atom command");

  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;

  if (atom->dpd_flag != 1)
    error->all(FLERR,
               "compute dpd requires atom_style with internal temperature and energies (e.g. dpd)");
}

/* ---------------------------------------------------------------------- */

ComputeDpdAtom::~ComputeDpdAtom()
{
  memory->destroy(dpdAtom);
}

/* ---------------------------------------------------------------------- */

void ComputeDpdAtom::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("^dpd/atom").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}

/* ----------------------------------------------------------------------
   gather compute vector data from other nodes
------------------------------------------------------------------------- */

void ComputeDpdAtom::compute_peratom()
{

  invoked_peratom = update->ntimestep;

  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *uChem = atom->uChem;
  double *dpdTheta = atom->dpdTheta;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  if (nlocal > nmax) {
    memory->destroy(dpdAtom);
    nmax = atom->nmax;
    memory->create(dpdAtom, nmax, size_peratom_cols, "dpd/atom:dpdAtom");
    array_atom = dpdAtom;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dpdAtom[i][0] = uCond[i];
      dpdAtom[i][1] = uMech[i];
      dpdAtom[i][2] = uChem[i];
      dpdAtom[i][3] = dpdTheta[i];
    } else {
      dpdAtom[i][0] = dpdAtom[i][1] = dpdAtom[i][2] = dpdAtom[i][3] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDpdAtom::memory_usage()
{
  double bytes = (double) size_peratom_cols * nmax * sizeof(double);
  return bytes;
}
