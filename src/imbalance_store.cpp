/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "imbalance_store.h"

#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceStore::ImbalanceStore(LAMMPS *lmp) : Imbalance(lmp), name(nullptr) {}

/* -------------------------------------------------------------------- */

ImbalanceStore::~ImbalanceStore()
{
  delete[] name;
}

/* -------------------------------------------------------------------- */

int ImbalanceStore::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal balance weight command");
  name = utils::strdup(arg[0]);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::compute(double *weight)
{
  int flag, cols;
  int index = atom->find_custom(name, flag, cols);

  // property does not exist

  if (index < 0 || flag != 1 || cols)
    error->all(FLERR, "Balance weight store vector does not exist");

  double *prop = atom->dvector[index];
  const int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; ++i) prop[i] = weight[i];
}

/* -------------------------------------------------------------------- */

std::string ImbalanceStore::info()
{
  return fmt::format("  storing weight in atom property d_{}\n", name);
}
