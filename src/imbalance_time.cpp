/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "pointers.h"
#include "imbalance_time.h"
#include "atom.h"
#include "error.h"
#include "comm.h"
#include "force.h"
#include "timer.h"

using namespace LAMMPS_NS;

int ImbalanceTime::options(int narg, char **arg)
{
  Error *error = _lmp->error;
  Force *force = _lmp->force;

  if (narg < 1) error->all(FLERR,"Illegal balance weight command");
  _factor = force->numeric(FLERR,arg[0]);
  if (_factor < 0.0 || _factor > 1.0)
    error->all(FLERR,"Illegal balance weight command");
  return 1;
}

/* -------------------------------------------------------------------- */

void ImbalanceTime::compute(double *weight)
{
  const int nlocal = _lmp->atom->nlocal;
  const int nprocs = _lmp->comm->nprocs;
  MPI_Comm world = _lmp->world;
  Timer *timer = _lmp->timer;

  if (_factor > 0.0) {
    // compute the cost function of based on relevant timers
    if (timer->has_normal()) {

      double cost = -_last;
      cost += timer->get_wall(Timer::PAIR);
      cost += timer->get_wall(Timer::NEIGH);
      cost += timer->get_wall(Timer::BOND);
      cost += timer->get_wall(Timer::KSPACE);

      double allcost;
      MPI_Allreduce(&cost,&allcost,1,MPI_DOUBLE,MPI_SUM,world);

      if (allcost > 0.0) {
        const double scale = (1.0-_factor) + _factor*cost*nprocs/allcost;
        for (int i = 0; i < nlocal; ++i) weight[i] *= scale;
      }
      // record time up to this point
      _last += cost;
    }
  }
}

/* -------------------------------------------------------------------- */

void ImbalanceTime::info(FILE *fp)
{
  if (_factor > 0.0)
    fprintf(fp,"  time weight factor: %g\n",_factor);
}
