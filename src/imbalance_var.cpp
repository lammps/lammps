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

#include "imbalance_var.h"

#include "atom.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "variable.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceVar::ImbalanceVar(LAMMPS *lmp) : Imbalance(lmp), name(nullptr) {}

/* -------------------------------------------------------------------- */

ImbalanceVar::~ImbalanceVar()
{
  delete[] name;
}

/* -------------------------------------------------------------------- */

int ImbalanceVar::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal balance weight command");
  name = utils::strdup(arg[0]);
  init(0);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::init(int /*flag*/)
{
  id = input->variable->find(name);
  if (id < 0) {
    error->all(FLERR, "Variable name for balance weight does not exist");
  } else {
    if (input->variable->atomstyle(id) == 0)
      error->all(FLERR, "Variable for balance weight has invalid style");
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::compute(double *weight)
{
  const int all = group->find("all");
  if (all < 0) return;

  double *values;
  const int nlocal = atom->nlocal;
  memory->create(values, nlocal, "imbalance:values");

  input->variable->compute_atom(id, all, values, 1, 0);

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (values[i] <= 0.0) flag = 1;
  int flagall;
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall) error->one(FLERR, "Balance weight <= 0.0");

  for (int i = 0; i < nlocal; i++) weight[i] *= values[i];

  memory->destroy(values);
}

/* -------------------------------------------------------------------- */

std::string ImbalanceVar::info()
{
  return fmt::format("  weight variable: {}\n", name);
}
