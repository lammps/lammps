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

   Harmonic restoring potential on the spin modulus,

       U = 1/2 k (|S| - S0)^2 ,

   which supplies the longitudinal energy scale that inertial spin dynamics
   needs.  The spin modulus is a free coordinate under fix nve/tspin, fix
   nvt/tspin and fix langevin/tspin, so unless the magnetic potential itself
   resolves the magnitude of the local moment and restores it, this fix or an
   equivalent one has to be present, otherwise the modulus is unbounded.
------------------------------------------------------------------------- */

#include "fix_spring_tspin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "respa.h"
#include "tspin.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

FixSpringTSpin::FixSpringTSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), k(0.0), s0(0.0), espring(0.0), ilevel_respa(0), index_sm(-1)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix spring/tspin", error);
  if (narg > 5) error->all(FLERR, 5, "Unknown fix spring/tspin keyword: {}", arg[5]);

  if (!atom->sp_flag)
    error->all(FLERR, "Fix spring/tspin requires atom style spin");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;

  k = utils::numeric(FLERR, arg[3], false, lmp);
  if (k <= 0.0) error->all(FLERR, 3, "Fix spring/tspin spring constant must be > 0");
  s0 = utils::numeric(FLERR, arg[4], false, lmp);
  if (s0 < 0.0) error->all(FLERR, 4, "Fix spring/tspin equilibrium modulus must be >= 0");
}

/* ---------------------------------------------------------------------- */

int FixSpringTSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::init()
{
  int flag, cols;
  index_sm = atom->find_custom(TSPIN_SMASS, flag, cols);
  if ((index_sm >= 0) && ((flag != 1) || (cols != 0))) index_sm = -1;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::setup(int vflag)
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

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::post_force(int /*vflag*/)
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  double *s_mass = (index_sm >= 0) ? atom->dvector[index_sm] : nullptr;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  espring = 0.0;
  const double hbar = force->hplanck / MY_2PI;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double smag = sp[i][3];
    const bool dynamic = s_mass && (s_mass[i] > 0.0);
    if ((smag <= TSPIN_EPS) && !dynamic) continue;

    const double dmod = smag - s0;
    espring += 0.5 * k * dmod * dmod;

    // The radial direction is undefined at exactly zero modulus.  Keep the
    // potential energy continuous for a dynamical spin and let its velocity
    // carry it through zero; the force is well-defined again immediately.

    if (smag == 0.0) continue;

    // F = -dU/dS = -k (|S|-S0) s^,  in eV/muB.  fm holds the precession
    // field of the SPIN package, so the force is stored the way pair style
    // deepspin stores it, scaled by |S|/hbar, and fix nve/tspin converts it
    // back.  A purely radial term like this one leaves fm x S unchanged and
    // therefore does not disturb the fixed-modulus styles.

    const double pref = -k * dmod * smag / hbar;
    fm[i][0] += pref * sp[i][0];
    fm[i][1] += pref * sp[i][1];
    fm[i][2] += pref * sp[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpringTSpin::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy stored in the spin modulus springs
------------------------------------------------------------------------- */

double FixSpringTSpin::compute_scalar()
{
  double all;
  MPI_Allreduce(&espring, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}
