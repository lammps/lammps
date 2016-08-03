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
#include "imbalance_var.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;

int ImbalanceVar::options(int narg, char **arg)
{
  Error *error = _lmp->error;
  Force *force = _lmp->force;

  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  int len = strlen(arg[0])+1;
  _name = new char[len];
  memcpy(_name,arg[0],len);
  this->init();

  return 1;
}

/* -------------------------------------------------------------------- */
 
void ImbalanceVar::init()
{
  Error *error = _lmp->error;
  Variable *variable = _lmp->input->variable;

  if (_name) {
    _id = variable->find(_name);
    if (_id < 0) {
      error->all(FLERR,"Variable name for balance weight does not exist");
    } else {
      if (variable->atomstyle(_id) == 0)
        error->all(FLERR,"Variable for balance weight has invalid style");
    }
  }
}

/* -------------------------------------------------------------------- */
 
void ImbalanceVar::compute(double *weight)
{
  if (_id >= 0) {
    const int all = _lmp->group->find("all");
    const int nlocal = _lmp->atom->nlocal;

    double *val = new double[nlocal];
    _lmp->input->variable->compute_atom(_id,all,val,1,0);
    for (int i = 0; i < nlocal; ++i) weight[i] *= val[i];
    delete[] val;
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceVar::info(FILE *fp)
{
  if (_id >= 0)
    fprintf(fp,"  weight variable: %s\n",_name);
}
