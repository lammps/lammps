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


#include "pointers.h"
#include "imbalance_group.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;

int ImbalanceGroup::options(int narg, char **arg)
{
  Error *error = _lmp->error;
  Force *force = _lmp->force;
  Group *group = _lmp->group;

  if (narg < 3) error->all(FLERR,"Illegal balance weight command");

  _num = force->inumeric(FLERR,arg[0]);
  if (_num < 1) error->all(FLERR,"Illegal balance weight command");
  if (2*_num+1 > narg) error->all(FLERR,"Illegal balance weight command");

  _id = new int[_num];
  _factor = new double[_num];
  for (int i = 0; i < _num; ++i) {
    _id[i] = group->find(arg[2*i+1]);
    if (_id[i] < 0)
      error->all(FLERR,"Unknown group in balance weight command");
    _factor[i] = force->numeric(FLERR,arg[2*i+2]);
  }
  return 2*_num+1;
}

/* -------------------------------------------------------------------- */

void ImbalanceGroup::compute(double *weight)
{
  const int * const mask = _lmp->atom->mask;
  const int * const bitmask = _lmp->group->bitmask;
  const int nlocal = _lmp->atom->nlocal;

  if (_num == 0) return;

  for (int i = 0; i < nlocal; ++i) {
    const int imask = mask[i];
    double iweight = weight[i];
    for (int j = 0; j < _num; ++j) {
      if (imask & bitmask[_id[j]])
        iweight *= _factor[j];
    }
    weight[i] = iweight;
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceGroup::info(FILE *fp)
{
  if (_num > 0) {
    const char * const * const names = _lmp->group->names;

    fprintf(fp,"  group weights:");
    for (int i = 0; i < _num; ++i)
      fprintf(fp," %s=%g",names[_id[i]],_factor[i]);
    fputs("\n",fp);
  }
}
