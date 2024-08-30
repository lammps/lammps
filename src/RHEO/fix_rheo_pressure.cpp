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
   Contributing authors:
   Joel Clemmer (SNL), Thomas O'Connor (CMU)
----------------------------------------------------------------------- */

#include "fix_rheo_pressure.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
enum { NONE, LINEAR, CUBIC, TAITWATER, TAITGENERAL };

static constexpr double SEVENTH = 1.0 / 7.0;

/* ---------------------------------------------------------------------- */

FixRHEOPressure::FixRHEOPressure(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), c_cubic(nullptr), csq(nullptr), csqinv(nullptr), rho0(nullptr),
    rho0inv(nullptr), tpower(nullptr), pbackground(nullptr), pressure_style(nullptr),
    fix_rheo(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix command");

  comm_forward = 1;

  // Currently can only have one instance of fix rheo/pressure
  if (igroup != 0) error->all(FLERR, "fix rheo/pressure command requires group all");

  int i, nlo, nhi;
  int n = atom->ntypes;
  memory->create(pressure_style, n + 1, "rheo:pressure_style");
  memory->create(c_cubic, n + 1, "rheo:c_cubic");
  memory->create(tpower, n + 1, "rheo:tpower");
  memory->create(pbackground, n + 1, "rheo:pbackground");
  for (i = 1; i <= n; i++) pressure_style[i] = NONE;

  int iarg = 3;
  while (iarg < narg) {
    utils::bounds(FLERR, arg[iarg], 1, n, nlo, nhi, error);

    if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/pressure", error);

    if (strcmp(arg[iarg + 1], "linear") == 0) {
      for (i = nlo; i <= nhi; i++) pressure_style[i] = LINEAR;
    } else if (strcmp(arg[iarg + 1], "tait/water") == 0) {
      for (i = nlo; i <= nhi; i++) pressure_style[i] = TAITWATER;
    } else if (strcmp(arg[iarg + 1], "tait/general") == 0) {
      if (iarg + 3 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/pressure tait", error);

      double tpower_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      double pbackground_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 2;

      for (i = nlo; i <= nhi; i++) {
        pressure_style[i] = TAITGENERAL;
        tpower[i] = tpower_one;
        pbackground[i] = pbackground_one;
      }
    } else if (strcmp(arg[iarg + 1], "cubic") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/pressure cubic", error);

      double c_cubic_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 1;

      for (i = nlo; i <= nhi; i++) {
        pressure_style[i] = CUBIC;
        c_cubic[i] = c_cubic_one;
      }
    } else {
      error->all(FLERR, "Illegal fix command, {}", arg[iarg]);
    }
    iarg += 2;
  }

  for (i = 1; i <= n; i++)
    if (pressure_style[i] == NONE)
      error->all(FLERR, "Must specify pressure for atom type {} in fix/rheo/pressure", i);
}

/* ---------------------------------------------------------------------- */

FixRHEOPressure::~FixRHEOPressure()
{
  memory->destroy(pressure_style);
  memory->destroy(csqinv);
  memory->destroy(rho0inv);
  memory->destroy(c_cubic);
  memory->destroy(tpower);
  memory->destroy(pbackground);
}

/* ---------------------------------------------------------------------- */

int FixRHEOPressure::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOPressure::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/pressure");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  csq = fix_rheo->csq;
  rho0 = fix_rheo->rho0;

  int n = atom->ntypes;
  memory->create(csqinv, n + 1, "rheo:rho0inv");
  memory->create(rho0inv, n + 1, "rheo:rho0inv");
  for (int i = 0; i <= n; i++) {
    csqinv[i] = 1.0 / csq[i];
    rho0inv[i] = 1.0 / rho0[i];
  }

  if (modify->get_fix_by_style("rheo/pressure").size() > 1)
    error->all(FLERR, "Can only specify one instance of fix rheo/pressure");
}

/* ---------------------------------------------------------------------- */

void FixRHEOPressure::setup_pre_force(int /*vflag*/)
{
  fix_rheo->pressure_fix_defined = 1;
  pre_force(0);
}

/* ----------------------------------------------------------------------
  Update (and forward) pressure every timestep
------------------------------------------------------------------------- */

void FixRHEOPressure::pre_force(int /*vflag*/)
{
  int *mask = atom->mask;
  int *type = atom->type;
  double *rho = atom->rho;
  double *pressure = atom->pressure;

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) pressure[i] = calc_pressure(rho[i], type[i]);

  if (comm_forward) comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int FixRHEOPressure::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                       int * /*pbc*/)
{
  double *pressure = atom->pressure;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = pressure[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOPressure::unpack_forward_comm(int n, int first, double *buf)
{
  double *pressure = atom->pressure;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) { pressure[i] = buf[m++]; }
}

/* ---------------------------------------------------------------------- */

double FixRHEOPressure::calc_pressure(double rho, int type)
{
  double p = 0.0;
  double dr, rr3, rho_ratio;

  if (pressure_style[type] == LINEAR) {
    p = csq[type] * (rho - rho0[type]);
  } else if (pressure_style[type] == CUBIC) {
    dr = rho - rho0[type];
    p = csq[type] * (dr + c_cubic[type] * dr * dr * dr);
  } else if (pressure_style[type] == TAITWATER) {
    rho_ratio = rho * rho0inv[type];
    rr3 = rho_ratio * rho_ratio * rho_ratio;
    p = csq[type] * rho0[type] * SEVENTH * (rr3 * rr3 * rho_ratio - 1.0);
  } else if (pressure_style[type] == TAITGENERAL) {
    rho_ratio = rho * rho0inv[type];
    p = csq[type] * rho0[type] * (pow(rho_ratio, tpower[type]) - 1.0) / tpower[type];
    p += pbackground[type];
  }
  return p;
}

/* ---------------------------------------------------------------------- */

double FixRHEOPressure::calc_rho(double p, int type)
{
  double rho = 0.0;

  if (pressure_style[type] == LINEAR) {
    rho = csqinv[type] * p + rho0[type];
  } else if (pressure_style[type] == CUBIC) {
    error->one(FLERR,
               "Rho calculation from pressure not yet supported for cubic pressure equation");
  } else if (pressure_style[type] == TAITWATER) {
    rho = pow(7.0 * p + csq[type] * rho0[type], SEVENTH);
    rho *= pow(rho0[type], 6.0 * SEVENTH);
    rho *= pow(csq[type], -SEVENTH);
  } else if (pressure_style[type] == TAITGENERAL) {
    p -= pbackground[type];
    rho = pow(tpower[type] * p + csq[type] * rho0[type], 1.0 / tpower[type]);
    rho *= pow(rho0[type], 1.0 - 1.0 / tpower[type]);
    rho *= pow(csq[type], -1.0 / tpower[type]);
  }
  return rho;
}
