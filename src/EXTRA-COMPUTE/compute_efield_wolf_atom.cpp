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

#include "compute_efield_wolf_atom.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PIS;

/* ---------------------------------------------------------------------- */

ComputeEfieldWolfAtom::ComputeEfieldWolfAtom(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), group2(nullptr), efield(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute efield/atom/wolf", error);

  peratom_flag = 1;
  size_peratom_cols = 3;

  nmax = -1;
  group2 = utils::strdup("all");
  jgroupbit = group->bitmask[0];
  cutoff_flag = 0;
  cutoff = 0.0;

  alpha = utils::numeric(FLERR, arg[3], false, lmp);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "group") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute efield/atom/wolf group", error);
      delete[] group2;
      group2 = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute efield/atom/wolf cutoff", error);
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      cutoff_flag = 1;
      iarg += 2;
    } else
      error->all(FLERR, "Unknown compute {} keyword: {}", style, arg[iarg]);
  }

  // sanity checks

  if (alpha <= 0.0) error->all(FLERR, "Compute efield/atom/wolf alpha value {} is invalid", alpha);
  if (cutoff_flag && cutoff <= 0.0)
    error->all(FLERR, "Compute efield/atom/wolf cutoff {} is invalid", cutoff);

  jgroup = group->find(group2);
  if (jgroup < 0) error->all(FLERR, "Compute efield/atom/wolf group {} does not exist", group2);
}

/* ---------------------------------------------------------------------- */

ComputeEfieldWolfAtom::~ComputeEfieldWolfAtom()
{
  delete[] group2;
  memory->destroy(efield);
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldWolfAtom::init()
{
  if (!atom->q_flag) error->all(FLERR, "Compute efield/wolf/atom requires atom attribute q");
  if (atom->mu_flag && (comm->me == 0))
    error->warning(FLERR, "Compute efield/wolf/atom does not support per-atom dipoles");

  // need an occasional full neighbor list
  auto req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  if (cutoff_flag) req->set_cutoff(cutoff);

  jgroup = group->find(group2);
  if (jgroup < 0) error->all(FLERR, "Compute efield/atom/wolf group {} does not exist", group2);
  jgroupbit = group->bitmask[jgroup];
}

/* ---------------------------------------------------------------------- */

void ComputeEfieldWolfAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

// clang-format off
/* ---------------------------------------------------------------------- */

void ComputeEfieldWolfAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow result array if necessary and clear

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(efield,nmax,3,"efield/atom/wolf:efield");
    array_atom = efield;
  }
  memset(&efield[0][0], 0, sizeof(double)*nmax*3);

  // invoke neighbor list build (will copy or build if necessary)

  neighbor->build_one(list);

  const auto inum = list->inum;
  const auto ilist = list->ilist;
  const auto numneigh = list->numneigh;
  const auto firstneigh = list->firstneigh;

  // compute coulomb force according to Wolf sum approximation
  const double * const * const x = atom->x;
  const int * const mask = atom->mask;
  const double * const q = atom->q;
  const double * const special_coul = force->special_coul;

  const double qqrd2e = force->qqrd2e;

  if (!cutoff_flag && force->pair) cutoff = force->pair->cutforce;
  const double cutsq = cutoff*cutoff;
  const double e_shift = erfc(alpha * cutoff) / cutoff;
  const double f_shift = -(e_shift + 2.0 * alpha / MY_PIS * exp(-alpha * alpha * cutsq)) / cutoff;

#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const auto jlist = firstneigh[i];
      const auto jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        const double factor_coul = special_coul[sbmask(j)];
        j &= NEIGHMASK;
        if (mask[j] & jgroupbit) {
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx*delx + dely*dely + delz*delz;
          if ((rsq > 0.0) && (rsq < cutsq)) {
            const double r = sqrt(rsq);
            double prefactor = qqrd2e * q[j] / r;
            double erfcc = erfc(alpha * r);
            double erfcd = exp(-alpha * alpha * r * r);
            double dvdrr = (erfcc / rsq + 2.0 * alpha / MY_PIS * erfcd / r) + f_shift;
            double forcecoul = dvdrr * rsq * prefactor;
            if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
            forcecoul /= rsq;

            efield[i][0] += delx * forcecoul;
            efield[i][1] += dely * forcecoul;
            efield[i][2] += delz * forcecoul;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEfieldWolfAtom::memory_usage()
{
  double bytes = 3.0 * nmax * sizeof(double);
  return bytes;
}
