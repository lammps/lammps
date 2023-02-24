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

#include "fix_alchemy.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

static double get_lambda(const bigint &step, const bigint &begin, const bigint &end, int iworld)
{
  double lambda = step - begin;
  if (lambda != 0.0) lambda /= end - begin;
  if (iworld == 0) lambda = 1.0 - lambda;
  return lambda;
}

/* ---------------------------------------------------------------------- */

FixAlchemy::FixAlchemy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), coordbuf(nullptr)
{
  lambda = epot[0] = epot[1] = 0.0;

  no_change_box = 1;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  // set up communicators
  if (universe->nworlds != 2) error->all(FLERR, "Must use exactly two partitions");
  int color = comm->me;
  int key = universe->iworld;
  MPI_Comm_split(universe->uworld, color, key, &samerank);

  if (narg != 3) error->all(FLERR, "Incorrect number of arguments for fix alchemy");
}

/* ---------------------------------------------------------------------- */

FixAlchemy::~FixAlchemy()
{
  MPI_Comm_free(&samerank);
}

/* ---------------------------------------------------------------------- */

int FixAlchemy::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::init() {}

/* ---------------------------------------------------------------------- */

void FixAlchemy::setup(int /*vflag*/)
{
  const int nall = atom->nlocal;
  const auto tag = atom->tag;
  const auto x = atom->x;

  memory->destroy(coordbuf);
  memory->create(coordbuf, 4 * sizeof(double) * nall, "alchemy:coordbuf");

  if (universe->iworld == 0) {
    int m = 0;
    for (int i = 0; i < nall; ++i) {
      coordbuf[m++] = ubuf(tag[i]).d;
      coordbuf[m++] = x[i][0];
      coordbuf[m++] = x[i][1];
      coordbuf[m++] = x[i][2];
    }
  }

  MPI_Bcast(coordbuf, 4 * nall, MPI_DOUBLE, 0, samerank);

  if (universe->iworld == 1) {
    int m = 0;
    for (int i = 0; i < nall; ++i) {
      tagint mytag = (tagint) ubuf(coordbuf[m++]).i;
      x[atom->map(mytag)][0] = coordbuf[m++];
      x[atom->map(mytag)][1] = coordbuf[m++];
      x[atom->map(mytag)][2] = coordbuf[m++];
    }
  }

  lambda = get_lambda(update->ntimestep, update->beginstep, update->endstep, universe->iworld);
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_integrate()
{
  int nall = 3 * (atom->nlocal + atom->nghost);
  //  MPI_Bcast(&atom->x[0][0], nall, MPI_DOUBLE, 0, samerank);
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_force(int /*vflag*/)
{
  lambda = get_lambda(update->ntimestep, update->beginstep, update->endstep, universe->iworld);
}

/* ---------------------------------------------------------------------- */

double FixAlchemy::compute_vector(int n)
{
  if (n == 0) return lambda;
  return epot[n - 1];
}
