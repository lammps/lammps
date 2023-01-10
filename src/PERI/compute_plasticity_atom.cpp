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

#include "compute_plasticity_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_peri_neigh.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePlasticityAtom::ComputePlasticityAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute plasticity/atom command");

  if (!force->pair_match("peri/eps",1))
    error->all(FLERR,"Compute plasticity/atom cannot be used "
               "with this pair style");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  plasticity = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputePlasticityAtom::~ComputePlasticityAtom()
{
  memory->destroy(plasticity);
}

/* ---------------------------------------------------------------------- */

void ComputePlasticityAtom::init()
{
  if ((comm->me == 0) && (modify->get_compute_by_style("plasticity/atom").size() > 1))
    error->warning(FLERR,"More than one compute plasticity/atom");

  // find associated PERI_NEIGH fix that must exist

  auto fixes = modify->get_fix_by_style("PERI_NEIGH");
  if (fixes.size() == 0)
    error->all(FLERR,"Compute plasticity/atom requires a peridynamic potential");
  else fix_peri_neigh = dynamic_cast<FixPeriNeigh *>(fixes.front());
}

/* ---------------------------------------------------------------------- */

void ComputePlasticityAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow damage array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(plasticity);
    nmax = atom->nmax;
    memory->create(plasticity,nmax,"plasticity/atom:plasticity");
    vector_atom = plasticity;
  }

  // extract plasticity for each atom in group

  double *lambdaValue = fix_peri_neigh->lambdaValue;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) plasticity[i] = lambdaValue[i];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePlasticityAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
