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

#include <string.h>
#include "pointers.h"
#include "imbalance_store.h"
#include "atom.h"
#include "error.h"
#include "input.h"

using namespace LAMMPS_NS;

int ImbalanceStore::options(int narg, char **arg)
{
  Error *error = _lmp->error;

  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  int len = strlen(arg[0])+1;
  _name = new char[len];
  memcpy(_name,arg[0],len);

  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::compute(double *weight)
{
  if (_name) {
    int dflag = 0;
    int idx = _lmp->atom->find_custom(_name,dflag);

    // property does not exist
    if (idx < 0 || dflag != 1) return;

    double *prop = _lmp->atom->dvector[idx];
    const int nlocal = _lmp->atom->nlocal;

    for (int i = 0; i < nlocal; ++i)
      prop[i] = weight[i];
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceStore::info(FILE *fp)
{
  if (_name)
    fprintf(fp,"  storing weight in atom property d_%s\n",_name);
}
