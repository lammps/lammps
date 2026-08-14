/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aram Davtyan (Deep Origin Inc.)
------------------------------------------------------------------------- */

#include "fix_baoab.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBAOAB::FixBAOAB(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), random(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix baoab", error);

  // --- Required arguments ---

  t_start = utils::numeric(FLERR, arg[3], false, lmp);
  t_stop = utils::numeric(FLERR, arg[4], false, lmp);
  damp = utils::numeric(FLERR, arg[5], false, lmp);
  seed = utils::inumeric(FLERR, arg[6], false, lmp);

  if (t_start < 0.0 || t_stop < 0.0) error->all(FLERR, 3, "Fix baoab temperatures must be >= 0");
  if (damp <= 0.0) error->all(FLERR, 5, "Fix baoab damp must be > 0");
  if (seed <= 0) error->all(FLERR, 6, "Fix baoab seed must be positive integer");

  // gamma = 1/damp (damp is the relaxation time, like fix langevin)
  gamma = 1.0 / damp;

  // --- Optional keyword arguments ---

  zeroflag = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "zero") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix baoab zero", error);
      zeroflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix baoab keyword {}", arg[iarg]);
    }
  }

  // --- Fix properties ---

  time_integrate = 1;         // this fix does time integration
  dynamic_group_allow = 1;    // works with dynamic groups
  scalar_flag = 1;            // produces a global scalar (energy)
  global_freq = 1;            // scalar computed every step
  extscalar = 1;              // scalar is extensive
  ecouple_flag = 1;           // outputs cumulative reservoir energy via compute_scalar()
  restart_global = 1;         // save/restore cumulative energy on restart

  // Initialize energy accounting
  energy = 0.0;
  energy_onestep = 0.0;

  // Set initial target
  t_target = t_start;

  // Each MPI rank gets a unique RNG stream
  random = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixBAOAB::~FixBAOAB()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixBAOAB::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;    // B A O A  (before force evaluation)
  mask |= FINAL_INTEGRATE;      // B        (after force evaluation)
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::init()
{
  // fix baoab and fix shake/rattle cannot be used together
  for (const auto &ifix : modify->get_fix_list())
    if ((utils::strmatch(ifix->style, "^shake") || utils::strmatch(ifix->style, "^rattle")) &&
        (ifix->groupbit & groupbit))
      error->all(FLERR, Error::NOLASTLINE, "Fix baoab is not compatible with fix {} on group {}",
                 ifix->style, group->names[igroup]);

  recompute_dt_coeffs();
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::recompute_dt_coeffs()
{
  double dt = update->dt;

  dtf = 0.5 * dt * force->ftm2v;    // half-kick: F/m -> dv
  dtby2 = 0.5 * dt;                 // half-drift: v -> dx
  c1 = exp(-gamma * dt);            // O-step decay (full dt)
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  double total = update->endstep - update->beginstep;
  if (total > 0)
    delta /= total;
  else
    delta = 0.0;

  t_target = t_start + delta * (t_stop - t_start);
}

/* ---------------------------------------------------------------------- */

// Called BEFORE force evaluation each timestep.
// Implements:  B (half kick) -> A (half drift) -> O (thermostat) -> A (half drift)
//
// The previous step's forces are still in f[] from final_integrate of
// the last step (or from setup on step 0).

void FixBAOAB::initial_integrate(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;    // per-atom mass (may be nullptr)
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // firstgroup optimization
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // Update target temperature (supports ramping via run start/stop)
  compute_target();

  // kT in velocity-squared * mass units
  //   force->boltz is kB in energy units of the chosen unit system
  //   force->mvv2e converts (1/2)mv^2 to energy: E = mvv2e * m * v^2
  //   so v^2 ~ kT / (mvv2e * m), hence noise scale = sqrt(kT / (mvv2e * m))
  double kT = force->boltz * t_target / force->mvv2e;

  // O-step coefficients that are mass-independent
  double one_minus_c1sq = 1.0 - c1 * c1;

  // For energy tally: accumulate KE change during O step
  energy_onestep = 0.0;

  // Accumulators for zeroing total random momentum
  double fran_total[3] = {0.0, 0.0, 0.0};
  double mass_total = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    // Handle both per-type and per-atom mass
    double mi = rmass ? rmass[i] : mass[type[i]];
    double invmass = 1.0 / mi;

    // ---- B step: half velocity kick from conservative forces ----
    //   dv = (F/m) * dt/2,  but in LAMMPS units:
    //   dv = ftm2v * (F/m) * dt/2 = dtf * F / m
    double dtfm = dtf * invmass;
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];

    // ---- A step: half position drift ----
    x[i][0] += dtby2 * v[i][0];
    x[i][1] += dtby2 * v[i][1];
    x[i][2] += dtby2 * v[i][2];

    // ---- O step: exact Ornstein-Uhlenbeck thermostat (full dt) ----
    //
    // Exact solution of  dv = -gamma*v*dt + sqrt(2*gamma*kT/m) dW:
    //   v_new = c1 * v_old + c2 * R,    R ~ N(0,1)
    //   c1 = exp(-gamma * dt)
    //   c2 = sqrt(kT / (mvv2e * m) * (1 - c1^2))
    //
    // Note: kT already divided by mvv2e above.

    double c2 = sqrt(kT * invmass * one_minus_c1sq);

    double ke_before =
        0.5 * force->mvv2e * mi * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

    double r0 = random->gaussian();
    double r1 = random->gaussian();
    double r2 = random->gaussian();

    v[i][0] = c1 * v[i][0] + c2 * r0;
    v[i][1] = c1 * v[i][1] + c2 * r1;
    v[i][2] = c1 * v[i][2] + c2 * r2;

    // Accumulate for zero-momentum correction
    if (zeroflag) {
      fran_total[0] += mi * c2 * r0;
      fran_total[1] += mi * c2 * r1;
      fran_total[2] += mi * c2 * r2;
      mass_total += mi;
    }

    double ke_after =
        0.5 * force->mvv2e * mi * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
    energy_onestep += ke_before - ke_after;

    // ---- A step: second half position drift ----
    x[i][0] += dtby2 * v[i][0];
    x[i][1] += dtby2 * v[i][1];
    x[i][2] += dtby2 * v[i][2];
  }

  // Zero net random momentum across all MPI ranks
  if (zeroflag) {
    double buf[4] = {fran_total[0], fran_total[1], fran_total[2], mass_total};
    double bufall[4];
    MPI_Allreduce(buf, bufall, 4, MPI_DOUBLE, MPI_SUM, world);

    if (bufall[3] > 0.0) {
      double inv_mtot = 1.0 / bufall[3];
      double vcorr[3];
      vcorr[0] = bufall[0] * inv_mtot;
      vcorr[1] = bufall[1] * inv_mtot;
      vcorr[2] = bufall[2] * inv_mtot;

      for (int i = 0; i < nlocal; i++) {
        if (!(mask[i] & groupbit)) continue;
        v[i][0] -= vcorr[0];
        v[i][1] -= vcorr[1];
        v[i][2] -= vcorr[2];

        // Also correct positions since the second A step already used
        // the un-corrected velocities
        x[i][0] -= dtby2 * vcorr[0];
        x[i][1] -= dtby2 * vcorr[1];
        x[i][2] -= dtby2 * vcorr[2];
      }
    }
  }

  energy += energy_onestep;
}

/* ---------------------------------------------------------------------- */

// Called AFTER force evaluation each timestep.
// Implements:  B (half kick with new forces)

void FixBAOAB::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // firstgroup optimization
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double mi = rmass ? rmass[i] : mass[type[i]];
    double dtfm = dtf / mi;

    // ---- B step: half velocity kick from new forces ----
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::reset_dt()
{
  recompute_dt_coeffs();
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::write_restart(FILE *fp)
{
  double energy_all;
  MPI_Allreduce(&energy, &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
  if (comm->me == 0) {
    int size = sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(&energy_all, sizeof(double), 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixBAOAB::restart(char *buf)
{
  // Rank 0 holds the global accumulated energy from before the restart;
  // other ranks start from 0. The Allreduce in compute_scalar sums them.
  energy = (comm->me == 0) ? *((double *) buf) : 0.0;
}

/* ---------------------------------------------------------------------- */

double FixBAOAB::compute_scalar()
{
  double energy_all;
  MPI_Allreduce(&energy, &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy_all;
}
