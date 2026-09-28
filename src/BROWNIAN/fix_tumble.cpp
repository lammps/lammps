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

#include "fix_tumble.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_const.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

enum { DIPOLE };

/* ---------------------------------------------------------------------- */

FixTumble::FixTumble(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), planar_flag(0), ptumble(0.0), nlevels_respa(1), rng(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix tumble", error);

  restart_global = 1;

  if (strcmp(arg[3], "dipole") == 0) {
    mode = DIPOLE;
  } else {
    error->all(FLERR, 3, "Unknown fix tumble mode {}", arg[3]);
  }

  rate = utils::numeric(FLERR, arg[4], false, lmp);
  if (rate <= 0.0) error->all(FLERR, 4, "Fix tumble rate must be > 0");

  seed = utils::inumeric(FLERR, arg[5], false, lmp);
  if (seed <= 0) error->all(FLERR, 5, "Fix tumble seed must be > 0");

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "planar_rotation") == 0) {
      if (domain->dimension == 2)
        error->all(FLERR, iarg,
                   "Fix tumble keyword planar_rotation is not allowed for 2d simulations");
      planar_flag = 1;
      ++iarg;
    } else {
      error->all(FLERR, iarg, "Unknown fix tumble keyword {}", arg[iarg]);
    }
  }

  // initialize Marsaglia RNG with processor-unique seed

  rng = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixTumble::~FixTumble()
{
  delete rng;
}

/* ---------------------------------------------------------------------- */

int FixTumble::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTumble::init()
{
  if ((mode == DIPOLE) && !atom->mu_flag)
    error->all(FLERR, Error::NOLASTLINE, "Fix tumble with mode dipole requires atom attribute mu");

  reset_dt();

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ----------------------------------------------------------------------
   tumbles are a Poisson process: probability of at least one tumble
   during one timestep for the given rate
------------------------------------------------------------------------- */

void FixTumble::reset_dt()
{
  ptumble = 1.0 - exp(-rate * update->dt);
}

/* ---------------------------------------------------------------------- */

void FixTumble::post_integrate()
{
  if (mode == DIPOLE) tumble_dipole();
}

/* ---------------------------------------------------------------------- */

void FixTumble::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_integrate();
}

/* ----------------------------------------------------------------------
   replace the dipole orientation of tumbling atoms with a random direction
   that is uniformly distributed on the unit circle (2d simulations or with
   planar_rotation) or on the unit sphere.  the dipole length is preserved.
------------------------------------------------------------------------- */

void FixTumble::tumble_dipole()
{
  double **mu = atom->mu;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int planar = (domain->dimension == 2) || planar_flag;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (rng->uniform() >= ptumble) continue;

    const double mulen = mu[i][3];
    if (mulen <= 0.0)
      error->one(FLERR, Error::NOLASTLINE,
                 "Fix tumble requires a non-zero dipole moment for atom {}", atom->tag[i]);

    const double phi = MY_2PI * rng->uniform();
    if (planar) {
      mu[i][0] = mulen * cos(phi);
      mu[i][1] = mulen * sin(phi);
      mu[i][2] = 0.0;
    } else {
      const double cz = 2.0 * rng->uniform() - 1.0;
      const double sz = sqrt(1.0 - cz * cz);
      mu[i][0] = mulen * sz * cos(phi);
      mu[i][1] = mulen * sz * sin(phi);
      mu[i][2] = mulen * cz;
    }
  }
}

/* ----------------------------------------------------------------------
   pack the per-processor RNG state into the restart file so that a run
   continued from a restart reproduces the original stochastic trajectory
------------------------------------------------------------------------- */

void FixTumble::write_restart(FILE *fp)
{
  int nsize = RanMars::STATE_SIZE * comm->nprocs + 1;    // pRNG state per proc + nprocs
  auto *list = new double[nsize];

  if (comm->me == 0) list[0] = comm->nprocs;

  double state[RanMars::STATE_SIZE];
  rng->get_state(state);
  MPI_Gather(state, RanMars::STATE_SIZE, MPI_DOUBLE, list + 1, RanMars::STATE_SIZE, MPI_DOUBLE, 0,
             world);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list, sizeof(double), nsize, fp);
  }
  delete[] list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restore the RNG state
------------------------------------------------------------------------- */

void FixTumble::restart(char *buf)
{
  auto *list = (double *) buf;

  int nprocs = (int) list[0];
  if (nprocs != comm->nprocs) {
    if (comm->me == 0)
      error->warning(FLERR, "Different number of procs. Cannot restore RNG state.");
  } else {
    // the size of the stored states depends on the version that wrote the restart file
    const int stride = RanMars::state_size(list + 1);
    rng->set_state(list + 1 + comm->me * stride);
  }
}
