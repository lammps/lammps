// clang-format off
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

#include "fix_heat_flow_sphere_temp.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {NONE, CONSTANT, TYPE};

/* ---------------------------------------------------------------------- */

FixHeatFlowSphereTemp::FixHeatFlowSphereTemp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  cp_style = NONE;

  int ntypes = atom->ntypes;
  if (strcmp(arg[3],"constant") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal fix command");
    cp_style = CONSTANT;
    cp = utils::numeric(FLERR,arg[4],false,lmp);
    if (cp < 0.0) error->all(FLERR,"Illegal fix command");
  } else if (strcmp(arg[3],"type") == 0) {
    if (narg != 4 + ntypes) error->all(FLERR,"Illegal fix command");
    cp_style = TYPE;
    memory->create(cp_type,ntypes+1,"fix/temp/integrate:cp_type");
    for (int i = 1; i <= ntypes; i++) {
      cp_type[i] = utils::numeric(FLERR,arg[3+i],false,lmp);
      if (cp_type[i] < 0.0) error->all(FLERR,"Illegal fix command");
    }
  } else {
    error->all(FLERR,"Illegal fix command");
  }

  if (cp_style == NONE)
    error->all(FLERR, "Must specify specific heat in fix temp/integrate");
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

int FixHeatFlowSphereTemp::setmask()
{
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlowSphereTemp::init()
{
  dt = update->dt;

  if (!atom->temperature_flag)
    error->all(FLERR,"Fix temp/integrate requires atom style with temperature property");
  if (!atom->heatflow_flag)
    error->all(FLERR,"Fix temp/integrate requires atom style with heatflow property");
}

/* ---------------------------------------------------------------------- */

void FixHeatFlowSphereTemp::final_integrate()
{
  // update temperature of atoms in group

  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        temperature[i] += dt * heatflow[i] / (calc_cp(i) * rmass[i]);
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        temperature[i] += dt * heatflow[i] / (calc_cp(i) * mass[type[i]]);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatFlowSphereTemp::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dt = update->dt;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixHeatFlowSphereTemp::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double FixHeatFlowSphereTemp::calc_cp(int i)
{
  if (cp_style == TYPE) {
    return cp_type[atom->type[i]];
  } else {
    return cp;
  }
}
