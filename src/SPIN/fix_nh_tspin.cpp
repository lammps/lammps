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

/* ------------------------------------------------------------------------
   Contributing author: AUTHOR_NAME_TBD (AFFILIATION_TBD)

   Nose-Hoover chain machinery for inertial spin dynamics.  This is a base
   class and registers no fix style of its own; fix nvt/tspin derives from
   it.  It does the same time integration as fix nve/tspin and in addition
   couples the spin velocities to a Nose-Hoover chain held at the requested
   spin temperature.  Only the spin degrees of freedom are thermostatted;
   use a separate fix (fix nvt, fix langevin, ...) if the lattice also has
   to be thermostatted.
------------------------------------------------------------------------- */

#include "fix_nh_tspin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNHTSpin::FixNHTSpin(LAMMPS *lmp, int narg, char **arg) :
    FixNVETSpin(lmp, narg, arg), tstat_flag(0), t_start(0.0), t_stop(0.0), t_target(0.0), t_freq(0.0),
    t_period(0.0), t_current(0.0), ke_target(0.0), tdof(0.0), drag(0.0), tdrag_factor(1.0),
    mtchain(3), nc_tchain(1), eta(nullptr), eta_dot(nullptr), eta_dotdot(nullptr),
    eta_mass(nullptr)
{
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  restart_global = 1;
  ecouple_flag = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temp") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " temp", error);
      t_start = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      t_stop = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      t_period = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (t_start < 0.0)
        error->all(FLERR, iarg + 1, "Fix {} start temperature must be >= 0", style);
      if (t_stop < 0.0) error->all(FLERR, iarg + 2, "Fix {} stop temperature must be >= 0", style);
      if (t_period <= 0.0)
        error->all(FLERR, iarg + 3, "Fix {} damping period must be > 0", style);
      tstat_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg], "tchain") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " tchain", error);
      mtchain = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (mtchain < 1) error->all(FLERR, iarg + 1, "Fix {} tchain must be >= 1", style);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tloop") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " tloop", error);
      nc_tchain = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (nc_tchain < 1) error->all(FLERR, iarg + 1, "Fix {} tloop must be >= 1", style);
      iarg += 2;
    } else if (strcmp(arg[iarg], "drag") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, std::string("fix ") + style + " drag", error);
      drag = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (drag < 0.0) error->all(FLERR, iarg + 1, "Fix {} drag must be >= 0", style);
      iarg += 2;

      // the remaining keywords belong to fix nve/tspin and are parsed there

    } else if (strcmp(arg[iarg], "lattice") == 0) {
      iarg += 2;
    } else if (strcmp(arg[iarg], "spin") == 0) {
      iarg += 2;
    } else if (strcmp(arg[iarg], "spinmass") == 0) {
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix {} keyword: {}", style, arg[iarg]);
    }
  }

  if (!spin_flag) error->all(FLERR, "Fix {} cannot be used with spin frozen", style);

  t_freq = tstat_flag ? 1.0 / t_period : 0.0;
  t_target = t_start;

  // Nose-Hoover chain state; eta_dot has one extra element so that the
  // chain update can always read eta_dot[ich+1]

  memory->create(eta, mtchain, "nh/tspin:eta");
  memory->create(eta_dot, mtchain + 1, "nh/tspin:eta_dot");
  memory->create(eta_dotdot, mtchain, "nh/tspin:eta_dotdot");
  memory->create(eta_mass, mtchain, "nh/tspin:eta_mass");

  for (int ich = 0; ich < mtchain; ich++) eta[ich] = eta_dot[ich] = eta_dotdot[ich] = 0.0;
  eta_dot[mtchain] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNHTSpin::~FixNHTSpin()
{
  memory->destroy(eta);
  memory->destroy(eta_dot);
  memory->destroy(eta_dotdot);
  memory->destroy(eta_mass);
}

/* ---------------------------------------------------------------------- */

int FixNHTSpin::setmask()
{
  return FixNVETSpin::setmask();
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::init()
{
  FixNVETSpin::init();
  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::reset_dt()
{
  FixNVETSpin::reset_dt();

  dthalf = 0.5 * update->dt;
  dt4 = 0.25 * update->dt;
  dt8 = 0.125 * update->dt;
  dtq = 0.5 * update->dt;

  tdrag_factor = 1.0 - (update->dt * t_freq * drag / nc_tchain);
}

/* ----------------------------------------------------------------------
   number of thermostatted spin degrees of freedom, and the chain masses
------------------------------------------------------------------------- */

void FixNHTSpin::setup(int vflag)
{
  FixNVETSpin::setup(vflag);

  double **sp = atom->sp;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  int ndof = 0;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && (sp[i][3] > 0.0) && (s_mass[i] > 0.0)) ndof += 3;

  int ndof_all = 0;
  MPI_Allreduce(&ndof, &ndof_all, 1, MPI_INT, MPI_SUM, world);
  tdof = ndof_all;

  if (tdof == 0.0) error->all(FLERR, "Fix {} has no spin degrees of freedom to thermostat", style);

  compute_temp_target();

  eta_mass[0] = tdof * force->boltz * t_target / (t_freq * t_freq);
  for (int ich = 1; ich < mtchain; ich++)
    eta_mass[ich] = force->boltz * t_target / (t_freq * t_freq);
  for (int ich = 1; ich < mtchain; ich++)
    eta_dotdot[ich] =
        (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - force->boltz * t_target) /
        eta_mass[ich];
}

/* ----------------------------------------------------------------------
   target spin temperature, ramped over the course of the run
------------------------------------------------------------------------- */

void FixNHTSpin::compute_temp_target()
{
  double delta = update->ntimestep - update->beginstep;
  if ((delta != 0.0) && (update->beginstep != update->endstep))
    delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
  ke_target = tdof * force->boltz * t_target;
}

/* ----------------------------------------------------------------------
   instantaneous spin temperature from the spin kinetic energy
------------------------------------------------------------------------- */

double FixNHTSpin::compute_spin_temp()
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double ke = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;
    ke += s_mass[i] *
        (s_dot[i][0] * s_dot[i][0] + s_dot[i][1] * s_dot[i][1] + s_dot[i][2] * s_dot[i][2]);
  }

  double keall = 0.0;
  MPI_Allreduce(&ke, &keall, 1, MPI_DOUBLE, MPI_SUM, world);

  if (tdof == 0.0) return 0.0;
  return keall * force->mvv2e / (tdof * force->boltz);
}

/* ----------------------------------------------------------------------
   scale all thermostatted spin velocities by factor
------------------------------------------------------------------------- */

void FixNHTSpin::nh_vs_scale(double factor)
{
  double **sp = atom->sp;
  double **s_dot = atom->v_s;
  double *s_mass = atom->s_mass;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if ((sp[i][3] <= 0.0) || (s_mass[i] <= 0.0)) continue;
    s_dot[i][0] *= factor;
    s_dot[i][1] *= factor;
    s_dot[i][2] *= factor;
  }
}

/* ----------------------------------------------------------------------
   propagate the Nose-Hoover chain by dthalf
   same scheme as FixNH::nhc_temp_integrate, acting on the spin velocities
------------------------------------------------------------------------- */

void FixNHTSpin::nhc_spin_integrate()
{
  int ich;
  double expfac;
  double kecurrent = tdof * force->boltz * t_current;

  eta_mass[0] = tdof * force->boltz * t_target / (t_freq * t_freq);
  for (ich = 1; ich < mtchain; ich++)
    eta_mass[ich] = force->boltz * t_target / (t_freq * t_freq);

  if (eta_mass[0] > 0.0)
    eta_dotdot[0] = (kecurrent - ke_target) / eta_mass[0];
  else
    eta_dotdot[0] = 0.0;

  const double ncfac = 1.0 / nc_tchain;
  for (int iloop = 0; iloop < nc_tchain; iloop++) {

    for (ich = mtchain - 1; ich > 0; ich--) {
      expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
      eta_dot[ich] *= expfac;
      eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
      eta_dot[ich] *= tdrag_factor;
      eta_dot[ich] *= expfac;
    }

    expfac = exp(-ncfac * dt8 * eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
    eta_dot[0] *= tdrag_factor;
    eta_dot[0] *= expfac;

    const double factor_eta = exp(-ncfac * dthalf * eta_dot[0]);
    nh_vs_scale(factor_eta);

    // rescale the temperature to account for the velocity scaling

    t_current *= factor_eta * factor_eta;
    kecurrent = tdof * force->boltz * t_current;

    if (eta_mass[0] > 0.0)
      eta_dotdot[0] = (kecurrent - ke_target) / eta_mass[0];
    else
      eta_dotdot[0] = 0.0;

    for (ich = 0; ich < mtchain; ich++) eta[ich] += ncfac * dthalf * eta_dot[ich];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac * dt4;
    eta_dot[0] *= expfac;

    for (ich = 1; ich < mtchain; ich++) {
      expfac = exp(-ncfac * dt8 * eta_dot[ich + 1]);
      eta_dot[ich] *= expfac;
      eta_dotdot[ich] =
          (eta_mass[ich - 1] * eta_dot[ich - 1] * eta_dot[ich - 1] - force->boltz * t_target) /
          eta_mass[ich];
      eta_dot[ich] += eta_dotdot[ich] * ncfac * dt4;
      eta_dot[ich] *= expfac;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::initial_integrate(int vflag)
{
  compute_temp_target();
  t_current = compute_spin_temp();
  nhc_spin_integrate();

  FixNVETSpin::initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::final_integrate()
{
  FixNVETSpin::final_integrate();

  t_current = compute_spin_temp();
  nhc_spin_integrate();
}

/* ----------------------------------------------------------------------
   energy of the Nose-Hoover chain, for the conserved quantity
------------------------------------------------------------------------- */

double FixNHTSpin::compute_scalar()
{
  double energy = tdof * force->boltz * t_target * eta[0];
  for (int ich = 1; ich < mtchain; ich++) energy += force->boltz * t_target * eta[ich];
  for (int ich = 0; ich < mtchain; ich++)
    energy += 0.5 * eta_mass[ich] * eta_dot[ich] * eta_dot[ich] * force->mvv2e;
  return energy;
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::write_restart(FILE *fp)
{
  const int nsize = 2 * mtchain + 1;
  auto list = new double[nsize];

  int n = 0;
  list[n++] = mtchain;
  for (int ich = 0; ich < mtchain; ich++) list[n++] = eta[ich];
  for (int ich = 0; ich < mtchain; ich++) list[n++] = eta_dot[ich];

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), nsize, fp);
  }
  delete[] list;
}

/* ---------------------------------------------------------------------- */

void FixNHTSpin::restart(char *buf)
{
  auto list = (double *) buf;

  int n = 0;
  const int mtchain_restart = static_cast<int>(list[n++]);
  if (mtchain_restart != mtchain)
    error->all(FLERR, "Fix {} tchain differs from the restart file value", style);

  for (int ich = 0; ich < mtchain; ich++) eta[ich] = list[n++];
  for (int ich = 0; ich < mtchain; ich++) eta_dot[ich] = list[n++];
}
