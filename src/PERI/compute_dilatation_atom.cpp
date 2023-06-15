// clang-format off
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
   Contributing author: Rezwanur Rahman, John Foster (UTSA)
------------------------------------------------------------------------- */

#include "compute_dilatation_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDilatationAtom::
ComputeDilatationAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute Dilatation/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  dilatation = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeDilatationAtom::~ComputeDilatationAtom()
{
  memory->destroy(dilatation);
}

/* ---------------------------------------------------------------------- */

void ComputeDilatationAtom::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("dilatation/atom").size() > 1))
    error->warning(FLERR,"More than one compute dilatation/atom");

  // check for compatible pair style

  if ((force->pair_match("^peri",0) == nullptr) || force->pair_match("^peri/pmb",0))
    error->all(FLERR,"Compute dilatation/atom cannot be used with this pair style");
}

/* ---------------------------------------------------------------------- */

void ComputeDilatationAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow dilatation array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(dilatation);
    nmax = atom->nmax;
    memory->create(dilatation,nmax,"dilatation/atom:dilatation");
    vector_atom = dilatation;
  }

  // extract dilatation for each atom in group

  int tmp;
  auto anypair = force->pair_match("^peri",0);
  auto theta = (double *)anypair->extract("theta",tmp);

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) dilatation[i] = theta[i];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDilatationAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
