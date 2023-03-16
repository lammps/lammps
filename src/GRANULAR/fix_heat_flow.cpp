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

#include "fix_heat_flow.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum {NONE, CONSTANT, TYPE};

/* ---------------------------------------------------------------------- */

FixHeatFlow::FixHeatFlow(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  cp_style = NONE;
  comm_forward = 1;
  comm_reverse = 1;

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

int FixHeatFlow::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::init()
{
  dt = update->dt;

  if (!atom->temperature_flag)
    error->all(FLERR,"Fix temp/integrate requires atom style with temperature property");
  if (!atom->heatflow_flag)
    error->all(FLERR,"Fix temp/integrate requires atom style with heatflow property");
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::setup(int /*vflag*/)
{
  // Identify whether this is the first instance of fix heat/flow
  first_flag = 0;

  int i = 0;
  auto fixlist = modify->get_fix_by_style("heat/flow");
  for (const auto &ifix : fixlist) {
    if (strcmp(ifix->id, id) == 0) break;
    i++;
  }

  if (i == 0) first_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::setup_pre_force(int /*vflag*/)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::pre_force(int /*vflag*/)
{
  // send updated temperatures to ghosts if first instance of fix
  // then clear heatflow for next force calculation
  double *heatflow = atom->heatflow;
  if (first_flag) {
    comm->forward_comm(this);
    for (int i = 0; i < atom->nmax; i++) heatflow[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::final_integrate()
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

  // add ghost contributions to heatflow if first instance of fix
  if (first_flag)
    comm->reverse_comm(this);

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

void FixHeatFlow::final_integrate_respa(int /*ilevel*/, int /*iloop*/)
{
  dt = update->dt;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double FixHeatFlow::calc_cp(int i)
{
  if (cp_style == TYPE) {
    return cp_type[atom->type[i]];
  } else {
    return cp;
  }
}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  double *temperature = atom->temperature;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = temperature[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;

  double *temperature = atom->temperature;

  for (i = first; i < last; i++) temperature[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  double *heatflow = atom->heatflow;

  for (int i = first; i < last; i++) {
    buf[m++] = heatflow[i];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  double *heatflow = atom->heatflow;

  for (int i = 0; i < n; i++)
    heatflow[list[i]] += buf[m++];
}
