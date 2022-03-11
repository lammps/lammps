/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "imbalance_time.h"

#include "atom.h"
#include "error.h"
#include "timer.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* -------------------------------------------------------------------- */

ImbalanceTime::ImbalanceTime(LAMMPS *lmp) : Imbalance(lmp) {}

/* -------------------------------------------------------------------- */

int ImbalanceTime::options(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal balance weight command");
  factor = utils::numeric(FLERR, arg[0], false, lmp);
  if (factor <= 0.0) error->all(FLERR, "Illegal balance weight command");
  return 1;
}

/* ----------------------------------------------------------------------
   reset last and timers if necessary
------------------------------------------------------------------------- */

void ImbalanceTime::init(int flag)
{
  last = 0.0;

  // flag = 1 if called from FixBalance at start of run
  //   init Timer, so accumulated time not carried over from previous run
  // should NOT init Timer if called from Balance, it uses time from last run

  if (flag) timer->init();
}

/* -------------------------------------------------------------------- */

void ImbalanceTime::compute(double *weight)
{
  if (!timer->has_normal()) return;

  // cost = CPU time for relevant timers since last invocation
  // localwt = weight assigned to each owned atom
  // just return if no time yet tallied
  // we 0.1 seconds as a minimum time to avoid computation of bogus
  // load balancing weights due to limited timer resolution/precision

  double cost = -last;
  cost += timer->get_wall(Timer::PAIR);
  cost += timer->get_wall(Timer::NEIGH);
  cost += timer->get_wall(Timer::BOND);
  cost += timer->get_wall(Timer::KSPACE);
  cost += 0.1;

  double maxcost;
  MPI_Allreduce(&cost, &maxcost, 1, MPI_DOUBLE, MPI_MAX, world);
  if (maxcost <= 0.1) return;

  int nlocal = atom->nlocal;
  double localwt = 0.0;
  if (nlocal) localwt = cost / nlocal;

  if (nlocal && localwt <= 0.0) error->one(FLERR, "Balance weight <= 0.0");

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

  // record time up to this point

  last += cost;
}

/* -------------------------------------------------------------------- */

std::string ImbalanceTime::info()
{
  return fmt::format("  time weight factor: {}\n", factor);
}
