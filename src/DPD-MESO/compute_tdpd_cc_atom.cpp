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

#include "compute_tdpd_cc_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTDPDCCAtom::ComputeTDPDCCAtom(LAMMPS *lmp, int narg, char **arg) : Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR, 2, "Number of arguments for compute tdpd/cc/atom command != 4");
  if (atom->tdpd_flag != 1)
    error->all(FLERR, 2, "compute tdpd/cc/atom command requires atom_style tdpd");

  index = utils::inumeric(FLERR, arg[3], false, lmp);

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  cc_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeTDPDCCAtom::~ComputeTDPDCCAtom()
{
  memory->sfree(cc_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeTDPDCCAtom::init()
{

  if ((comm->me == 0) && (modify->get_compute_by_style("^tdpd/cc/atom").size() > 1))
    error->warning(FLERR, "More than one compute {}", style);
}

/* ---------------------------------------------------------------------- */

void ComputeTDPDCCAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow cc_vector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(cc_vector);
    nmax = atom->nmax;
    cc_vector = (double *) memory->smalloc(nmax * sizeof(double), "tdpd/cc/atom:cc_vector");
    vector_atom = cc_vector;
  }

  double **cc = atom->cc;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cc_vector[i] = cc[i][index - 1];
    } else
      cc_vector[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeTDPDCCAtom::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
