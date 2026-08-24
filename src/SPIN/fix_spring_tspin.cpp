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

   Harmonic restoring potential on the spin modulus,

       U = 1/2 k (|S| - S0)^2 ,

   which supplies the longitudinal energy scale that inertial spin dynamics
   needs.  The SPIN pair styles provide a Landau-Lifshitz effective field,
   whose component parallel to S does no work on a fixed-modulus spin but
   acts as a radial driving force once the modulus is dynamical.  A run
   using fix nve/tspin, fix nvt/tspin or fix langevin/tspin therefore has to
   include this fix (or another longitudinal potential), otherwise the spin
   modulus is unbounded.
------------------------------------------------------------------------- */

#include "fix_spring_tspin.h"

#include "atom.h"
#include "error.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpringTSpin::FixSpringTSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), k(0.0), s0(0.0), espring(0.0), ilevel_respa(0)
{
  if (narg != 5) utils::missing_cmd_args(FLERR, "fix spring/tspin", error);

  if (!atom->tsp_flag)
    error->all(FLERR, "Fix spring/tspin requires atom style tspin");

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
  double **f_spin = atom->f_spin;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  espring = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double smag = sp[i][3];
    if (smag <= 0.0) continue;

    const double dmod = smag - s0;
    espring += 0.5 * k * dmod * dmod;

    // F = -dU/dS = -k (|S|-S0) S/|S|, and S/|S| is the spin direction
    // this is a genuine energy gradient, so it goes into f_spin rather than
    // into the precession field fm

    const double pref = -k * dmod;
    f_spin[i][0] += pref * sp[i][0];
    f_spin[i][1] += pref * sp[i][1];
    f_spin[i][2] += pref * sp[i][2];
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
