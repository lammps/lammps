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

#include "fix_dt_reset.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "integrate.h"
#include "lattice.h"
#include "modify.h"
#include "output.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

FixDtReset::FixDtReset(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR, "Illegal fix dt/reset command");

  // set time_depend, else elapsed time accumulation can be messed up

  time_depend = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;
  dynamic_group_allow = 1;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix dt/reset command");

  minbound = maxbound = 1;
  tmin = tmax = 0.0;
  if (strcmp(arg[4], "NULL") == 0)
    minbound = 0;
  else
    tmin = utils::numeric(FLERR, arg[4], false, lmp);
  if (strcmp(arg[5], "NULL") == 0)
    maxbound = 0;
  else
    tmax = utils::numeric(FLERR, arg[5], false, lmp);
  xmax = utils::numeric(FLERR, arg[6], false, lmp);

  if (minbound && tmin < 0.0) error->all(FLERR, "Illegal fix dt/reset command");
  if (maxbound && tmax < 0.0) error->all(FLERR, "Illegal fix dt/reset command");
  if (minbound && maxbound && tmin >= tmax) error->all(FLERR, "Illegal fix dt/reset command");
  if (xmax <= 0.0) error->all(FLERR, "Illegal fix dt/reset command");

  int scaleflag = 1;

  emax = -1.0;
  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix dt/reset command");
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Illegal fix dt/reset command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "emax") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix dt/reset command");
      emax = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (emax <= 0.0) error->all(FLERR, "Illegal fix dt/reset command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix dt/reset command");
  }

  // setup scaling, based on xlattice parameter

  if (scaleflag) xmax *= domain->lattice->xlattice;

  // initializations

  t_laststep = 0.0;
  laststep = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixDtReset::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDtReset::init()
{
  // set rRESPA flag

  respaflag = 0;
  if (utils::strmatch(update->integrate_style, "^respa")) respaflag = 1;

  ftm2v = force->ftm2v;
  mvv2e = force->mvv2e;
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixDtReset::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixDtReset::end_of_step()
{
  double dtv, dtf, dte, dtsq;
  double vsq, fsq, massinv;
  double delx, dely, delz, delr;

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dtmin = BIG;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass)
        massinv = 1.0 / rmass[i];
      else
        massinv = 1.0 / mass[type[i]];
      vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
      fsq = f[i][0] * f[i][0] + f[i][1] * f[i][1] + f[i][2] * f[i][2];
      dtv = dtf = dte = BIG;
      if (vsq > 0.0) dtv = xmax / sqrt(vsq);
      if (fsq > 0.0) dtf = sqrt(2.0 * xmax / (ftm2v * sqrt(fsq) * massinv));
      dt = MIN(dtv, dtf);
      if ((emax > 0.0) && (fsq * vsq > 0.0)) {
        dte = emax / sqrt(fsq * vsq) / sqrt(ftm2v * mvv2e);
        dt = MIN(dt, dte);
      }
      dtsq = dt * dt;
      delx = dt * v[i][0] + 0.5 * dtsq * massinv * f[i][0] * ftm2v;
      dely = dt * v[i][1] + 0.5 * dtsq * massinv * f[i][1] * ftm2v;
      delz = dt * v[i][2] + 0.5 * dtsq * massinv * f[i][2] * ftm2v;
      delr = sqrt(delx * delx + dely * dely + delz * delz);
      if (delr > xmax) dt *= xmax / delr;
      dtmin = MIN(dtmin, dt);
    }

  MPI_Allreduce(&dtmin, &dt, 1, MPI_DOUBLE, MPI_MIN, world);

  if (minbound) dt = MAX(dt, tmin);
  if (maxbound) dt = MIN(dt, tmax);

  // if timestep didn't change, just return
  // else reset update->dt and other classes that depend on it
  // rRESPA, pair style, fixes

  if (dt == update->dt) return;

  laststep = update->ntimestep;

  // calls to other classes that need to know timestep size changed
  // similar logic is in Input::timestep()

  update->update_time();
  update->dt = dt;
  update->dt_default = 0;
  if (respaflag) update->integrate->reset_dt();
  if (force->pair) force->pair->reset_dt();
  for (auto &ifix : modify->get_fix_list()) ifix->reset_dt();
  output->reset_dt();
}

/* ---------------------------------------------------------------------- */

double FixDtReset::compute_scalar()
{
  return (double) laststep;
}
