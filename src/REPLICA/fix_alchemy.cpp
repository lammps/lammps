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
#include "compute.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
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

FixAlchemy::FixAlchemy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR, "Incorrect number of arguments for fix alchemy");

  lambda = epot[0] = epot[1] = 0.0;

  no_change_box = 1;
  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  ilevel_respa = 0;

  // set up communicators
  if (universe->nworlds != 2) error->all(FLERR, "Must use exactly two partitions");
  int color = comm->me;
  int key = universe->iworld;
  MPI_Comm_split(universe->uworld, color, key, &samerank);

  id_pe = std::string(id) + "_pe";
  modify->add_compute(id_pe + " all pe");
}

/* ---------------------------------------------------------------------- */

FixAlchemy::~FixAlchemy()
{
  MPI_Comm_free(&samerank);
  modify->delete_compute(id_pe);
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

void FixAlchemy::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    auto respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    respa->copy_f_flevel(ilevel_respa);
  } else {
    post_force(vflag);
  }
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_integrate()
{
  const int nall = atom->nlocal + atom->nghost;
  MPI_Bcast(&atom->x[0][0], 3 * nall, MPI_DOUBLE, 0, samerank);
}

/* ---------------------------------------------------------------------- */

void FixAlchemy::post_force(int /*vflag*/)
{
  const int nall = atom->nlocal;
  auto f = atom->f;
  lambda = get_lambda(update->ntimestep, update->beginstep, update->endstep, universe->iworld);

  double *commbuf = new double[3 * nall];
  int m = 0;
  for (int i = 0; i < nall; ++i) {
    commbuf[m++] = f[i][0] * lambda;
    commbuf[m++] = f[i][1] * lambda;
    commbuf[m++] = f[i][2] * lambda;
  }
  MPI_Allreduce(commbuf, &f[0][0], 3 * nall, MPI_DOUBLE, MPI_SUM, samerank);

  auto pe = modify->get_compute_by_id(id_pe);
  commbuf[0] = commbuf[1] = commbuf[2] = 0.0;
  commbuf[universe->iworld] = pe->compute_scalar() / comm->nprocs;
  commbuf[2] = lambda * pe->compute_scalar() / comm->nprocs;
  MPI_Allreduce(commbuf, epot, 3, MPI_DOUBLE, MPI_SUM, universe->uworld);
  delete[] commbuf;
  pe->addstep(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

double FixAlchemy::compute_vector(int n)
{
  if (n == 0) return lambda;
  return epot[n - 1];
}
