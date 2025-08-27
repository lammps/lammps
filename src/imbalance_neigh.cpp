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

#include "imbalance_neigh.h"

#include "accelerator_kokkos.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "neighbor.h"

using namespace LAMMPS_NS;

static constexpr double BIG = 1.0e20;

/* -------------------------------------------------------------------- */

ImbalanceNeigh::ImbalanceNeigh(LAMMPS *lmp) : Imbalance(lmp)
{
  did_warn = 0;
}

/* -------------------------------------------------------------------- */

int ImbalanceNeigh::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal balance weight command");
  factor = utils::numeric(FLERR, arg[0], false, lmp);
  if (factor <= 0.0) error->all(FLERR, "Illegal balance weight command");
  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceNeigh::compute(double *weight)
{
  if (factor == 0.0) return;

  // cannot use neighbor list weight with KOKKOS using GPUs

  if (lmp->kokkos && lmp->kokkos->kokkos_exists) {
    if (lmp->kokkos->ngpus > 0) {
      if (comm->me == 0 && !did_warn)
        error->warning(FLERR, "Balance weight neigh skipped with KOKKOS using GPUs");
      did_warn = 1;
      return;
    }
  }

  bigint neighsum = neighbor->get_nneigh_half();
  if (neighsum < 0) neighsum = neighbor->get_nneigh_full();

  if ((neighsum < 0) || (neighbor->ago < 0)) {
    if (comm->me == 0 && !did_warn)
      error->warning(FLERR, "Balance weight neigh skipped b/c no suitable list found");
    did_warn = 1;
    return;
  }

  // neighsum = total neigh count for atoms on this proc
  // localwt = weight assigned to each owned atom

  double localwt = 0.0;
  const int nlocal = atom->nlocal;
  if (nlocal) localwt = 1.0 * neighsum / nlocal;

  if (nlocal && localwt < 0.0) error->one(FLERR, "Balance weight < 0.0");

  // apply factor if specified != 1.0
  // wtlo,wthi = lo/hi values excluding 0.0 due to no atoms on this proc
  // lo value does not change
  // newhi = new hi value to give hi/lo ratio factor times larger/smaller
  // expand/contract all localwt values from lo->hi to lo->newhi

  if (factor != 1.0) {
    double wtlo, wthi;
    if (localwt == 0.0) localwt = BIG;
    MPI_Allreduce(&localwt, &wtlo, 1, MPI_DOUBLE, MPI_MIN, world);
    if (localwt == BIG) localwt = 0.0;
    MPI_Allreduce(&localwt, &wthi, 1, MPI_DOUBLE, MPI_MAX, world);
    if (wtlo == wthi) return;

    double newhi = wthi * factor;
    localwt = wtlo + ((localwt - wtlo) / (wthi - wtlo)) * (newhi - wtlo);
  }

  for (int i = 0; i < nlocal; i++) weight[i] *= localwt;
}

/* -------------------------------------------------------------------- */

std::string ImbalanceNeigh::info()
{
  return fmt::format("  neighbor weight factor: {}\n", factor);
}
