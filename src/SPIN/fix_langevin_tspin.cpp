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
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com

   Langevin thermostat for inertial spin dynamics.  The full bath force on
   every spin in the group is

       F_bath = -(s_mass/t_period) v_s
                + sqrt(2 s_mass k_B T/(t_period dt)) R

   and is encoded in the rad.THz fm array as fm += |S|/hbar * F_bath.

   so that the spin velocities v_s sample a Maxwell-Boltzmann distribution at
   the requested spin temperature.  This is a kinetic thermostat acting on
   the spin momenta, and is unrelated to fix langevin/spin, which adds
   transverse Gilbert damping and a fluctuating field to the fixed-modulus
   Landau-Lifshitz-Gilbert equation.
------------------------------------------------------------------------- */

#include "fix_langevin_tspin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "random_mars.h"
#include "respa.h"
#include "tspin.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

FixLangevinTSpin::FixLangevinTSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), t_target(0.0), zeroflag(0), index_vs(-1), index_sm(-1), ilevel_respa(0), gamma_drag(0.0),
    gamma_random(0.0), random(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix langevin/tspin", error);

  if (!atom->sp_flag)
    error->all(FLERR, "Fix langevin/tspin requires atom style spin");

  respa_level_support = 1;

  t_start = utils::numeric(FLERR, arg[3], false, lmp);
  t_stop = utils::numeric(FLERR, arg[4], false, lmp);
  if (t_start < 0.0) error->all(FLERR, 3, "Fix langevin/tspin start temperature must be >= 0");
  if (t_stop < 0.0) error->all(FLERR, 4, "Fix langevin/tspin stop temperature must be >= 0");
  t_target = t_start;

  t_period = utils::numeric(FLERR, arg[5], false, lmp);
  if (t_period <= 0.0) error->all(FLERR, 5, "Fix langevin/tspin damping period must be > 0");

  const int seed = utils::inumeric(FLERR, arg[6], false, lmp);
  if (seed <= 0) error->all(FLERR, 6, "Fix langevin/tspin seed must be > 0");

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "zero") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix langevin/tspin zero", error);
      zeroflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix langevin/tspin keyword: {}", arg[iarg]);
    }
  }

  random = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixLangevinTSpin::~FixLangevinTSpin()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixLangevinTSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinTSpin::init()
{
  tspin_bind_state(atom, error, std::string("Fix ") + style, index_vs, index_sm);
  reset_dt();

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevinTSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    auto respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    respa->copy_f_flevel(ilevel_respa);
  }
}

/* ----------------------------------------------------------------------
   target spin temperature, ramped over the course of the run
------------------------------------------------------------------------- */

void FixLangevinTSpin::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if ((delta != 0.0) && (update->beginstep != update->endstep))
    delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop - t_start);
}

/* ---------------------------------------------------------------------- */

void FixLangevinTSpin::post_force(int /*vflag*/)
{
  compute_target();

  double **sp = atom->sp;
  double **fm = atom->fm;
  double **s_dot = atom->darray[index_vs];
  double *s_mass = atom->dvector[index_sm];
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  const double tsqrt = sqrt(t_target);
  const double hbar = force->hplanck / MY_2PI;

  // the friction and noise prefactors are per atom, because the spin mass is
  // a per-atom property

  double fsum[3] = {0.0, 0.0, 0.0};
  int count = 0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (s_mass[i] <= 0.0) continue;

    const double gamma1 = -s_mass[i] * gamma_drag;
    const double gamma2 = sqrt(s_mass[i]) * gamma_random * tsqrt;

    double fran[3];
    fran[0] = gamma2 * (random->uniform() - 0.5);
    fran[1] = gamma2 * (random->uniform() - 0.5);
    fran[2] = gamma2 * (random->uniform() - 0.5);

    // the bath force is an energy gradient in eV/muB; fm keeps the units of
    // the SPIN package, so scale by |S|/hbar on the way in

    const double tofm = sp[i][3] / hbar;
    fm[i][0] += tofm * (gamma1 * s_dot[i][0] + fran[0]);
    fm[i][1] += tofm * (gamma1 * s_dot[i][1] + fran[1]);
    fm[i][2] += tofm * (gamma1 * s_dot[i][2] + fran[2]);

    if (zeroflag) {
      fsum[0] += fran[0];
      fsum[1] += fran[1];
      fsum[2] += fran[2];
      count++;
    }
  }

  // optionally remove the net random force on the spin system

  if (!zeroflag) return;

  double fsumall[3];
  int countall;
  MPI_Allreduce(fsum, fsumall, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, world);
  if (countall == 0) return;

  fsumall[0] /= countall;
  fsumall[1] /= countall;
  fsumall[2] /= countall;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (s_mass[i] <= 0.0) continue;
    const double tofm = sp[i][3] / hbar;
    fm[i][0] -= tofm * fsumall[0];
    fm[i][1] -= tofm * fsumall[1];
    fm[i][2] -= tofm * fsumall[2];
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevinTSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLangevinTSpin::reset_dt()
{
  // F = -(m_s/t_period) v_s + sqrt(24 kB T m_s/(t_period dt)) (U[0,1) - 1/2)
  // the uniform deviate has variance 1/12, hence the factor 24
  // the prefactors are those of the lattice fix langevin, so the force comes
  // out in the same units as an ordinary force and post_force() converts it

  gamma_drag = 1.0 / t_period / force->ftm2v;
  gamma_random = sqrt(24.0 * force->boltz / t_period / update->dt / force->mvv2e) / force->ftm2v;
}
