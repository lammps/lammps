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
#include "imbalance_store.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* -------------------------------------------------------------------- */

ImbalanceStore::ImbalanceStore(LAMMPS *lmp) : Imbalance(lmp), name(0) {}

/* -------------------------------------------------------------------- */

ImbalanceStore::~ImbalanceStore()
{
  delete [] name;
}

/* -------------------------------------------------------------------- */

int ImbalanceStore::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal balance weight command");

  int len = strlen(arg[0]) + 1;
  name = new char[len];
  memcpy(name,arg[0],len);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::compute(double *weight)
{
  int dflag = 0;
  int idx = atom->find_custom(name,dflag);

  // property does not exist

  if (idx < 0 || dflag != 1) return;

  double *prop = atom->dvector[idx];
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i)
    prop[i] = weight[i];
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::info(FILE *fp)
{
  fprintf(fp,"  storing weight in atom property d_%s\n",name);
}
