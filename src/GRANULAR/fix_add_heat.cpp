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
   Contributing author: Joel Clemmer (SNL)
------------------------------------------------------------------------- */

#include "fix_add_heat.h"

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { CONSTANT, EQUAL, ATOM };
enum { ADD, LINEAR, QUARTIC };

/* ---------------------------------------------------------------------- */

FixAddHeat::FixAddHeat(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), varstr(nullptr), vatom(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix add/heat", error);
  dynamic_group_allow = 1;
  overwrite_flag = 0;

  if (strcmp(arg[3], "constant") == 0) {
    style = ADD;
  } else if (strcmp(arg[3], "linear") == 0) {
    style = LINEAR;
  } else if (strcmp(arg[3], "quartic") == 0) {
    style = QUARTIC;
  } else {
    error->all(FLERR, "Invalid option {}", arg[3]);
  }

  if (utils::strmatch(arg[4], "^v_")) {
    varstr = utils::strdup(arg[4] + 2);
  } else {
    value = utils::numeric(FLERR, arg[4], false, lmp);
    vstyle = CONSTANT;
  }

  int iarg = 5;
  if (style != ADD) {
    if (narg != 6) utils::missing_cmd_args(FLERR, "fix add/heat", error);
    prefactor = utils::numeric(FLERR, arg[5], false, lmp);
    iarg = 6;
  }

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg], "overwrite") == 0) {
      if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix add/heat", error);
      overwrite_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix add/heat command, invalid argument {}", arg[iarg]);
    }
  }

  maxatom = -1;
}

/* ---------------------------------------------------------------------- */

FixAddHeat::~FixAddHeat()
{
  delete[] varstr;
  memory->destroy(vatom);
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
      vstyle = EQUAL;
    else if (input->variable->atomstyle(var))
      vstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix addforce is invalid style", varstr);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddHeat::post_force(int /*vflag*/)
{
  int *mask = atom->mask;
  double *heatflow = atom->heatflow;
  double *temperature = atom->temperature;

  if (vstyle == ATOM) {
    if (atom->nmax > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(vatom);
      memory->create(vatom, maxatom, "addheat:vatom");
    }

    input->variable->compute_atom(var, igroup, &vatom[0], 1, 0);
  }

  if (overwrite_flag)
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        heatflow[i] = 0.0;

  double vtmp = 0.0;
  if (vstyle == CONSTANT) vtmp = value;
  if (vstyle == EQUAL) vtmp = input->variable->compute_equal(var);
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (vstyle == ATOM) vtmp = vatom[i];

      if (style == ADD) {
        heatflow[i] += vtmp;
      } else if (style == LINEAR) {
        heatflow[i] += prefactor * (vtmp - temperature[i]);
      } else if (style == QUARTIC) {
        heatflow[i] += prefactor * (pow(vtmp, 4.0) - pow(temperature[i], 4.0));
      }
    }
  }
}
