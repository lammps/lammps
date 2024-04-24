/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_add_heat.h"

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixAddHeat::FixAddHeat(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), varstr(nullptr), qatom(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix add/heat", error);
  dynamic_group_allow = 1;
  overwrite_flag = 0;

  style = NONE;
  if (utils::strmatch(arg[3], "^v_")) {
    varstr = utils::strdup(arg[3] + 2);
  } else {
    value = utils::numeric(FLERR, arg[3], false, lmp);
    style = CONSTANT;
  }

  // optional args

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "overwrite") == 0) {
      overwrite_flag = 1;
      iarg += 1;
    } else
      error->all(FLERR, "Illegal fix viscous command");
  }

  maxatom = -1;
}

/* ---------------------------------------------------------------------- */

FixAddHeat::~FixAddHeat()
{
  delete[] varstr;
  memory->destroy(qatom);
}

/* ---------------------------------------------------------------------- */

int FixAddHeat::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddHeat::init()
{
  if (!atom->temperature_flag)
    error->all(FLERR, "Fix add/heat requires atom style with temperature property");
  if (!atom->heatflow_flag)
    error->all(FLERR, "Fix add/heat requires atom style with heatflow property");

  // check variable

  if (varstr) {
    var = input->variable->find(varstr);
    if (var < 0) error->all(FLERR, "Variable {} for fix addforce does not exist", varstr);
    if (input->variable->equalstyle(var))
      style = EQUAL;
    else if (input->variable->atomstyle(var))
      style = ATOM;
    else
      error->all(FLERR, "Variable {} for fix addforce is invalid style", varstr);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddHeat::post_force(int /*vflag*/)
{
  int *mask = atom->mask;
  double *heatflow = atom->heatflow;
  double dtinv = 1.0 / update->dt;

  if (overwrite_flag) {
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        heatflow[i] = 0.0;
  }

  if (style == CONSTANT) {
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        heatflow[i] += value * dtinv;
  } else if (style == EQUAL) {
    value = input->variable->compute_equal(var);
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        heatflow[i] += value * dtinv;
  } else if (style == ATOM) {

    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(qatom);
      memory->create(qatom, maxatom, "addheat:qatom");
    }
    input->variable->compute_atom(var, igroup, &qatom[0], 1, 0);
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        heatflow[i] += qatom[i] * dtinv;
  }
}
