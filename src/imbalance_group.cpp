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

int ImbalanceGroup::options(LAMMPS *lmp, int narg, char **arg)
{
  Error *error = lmp->error;
  Force *force = lmp->force;
  Group *group = lmp->group;

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
  return _num;
}
 
void ImbalanceGroup::compute(LAMMPS *lmp, double *weight)
{
  const int * const mask = lmp->atom->mask;
  const int * const bitmask = lmp->group->bitmask;
  const int nlocal = lmp->atom->nlocal;

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
