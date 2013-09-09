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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_countdown.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "update.h"

#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

static const double decay = 0.999;
static const int maxtimes = 1000;

/* ---------------------------------------------------------------------- */

FixCountdown::FixCountdown(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 4;

  _times = NULL;
  _last_wall = _complete = 0.0;
  _outfreq = 0;
  _ntimes = _idx = 0;

  if (narg == 5) {
    nevery = force->inumeric(FLERR,arg[3]);
    _outfreq = force->inumeric(FLERR,arg[4]);
  } else error->all(FLERR,"Illegal fix countdown command");

  if (nevery <= 0)
    error->all(FLERR,"Fix countdown requires: Nevery > 0");

  global_freq = nevery;
}

/* ---------------------------------------------------------------------- */

FixCountdown::~FixCountdown()
{
  memory->destroy(_times);
}

/* ---------------------------------------------------------------------- */

int FixCountdown::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCountdown::setup(int vflag)
{
  if (_times == NULL) {
    memory->create(_times,maxtimes,"countdown:times");
  }
  _ntimes = _idx = 0;
  _last_wall = MPI_Wtime();
  _complete = 0.0;
  _eta[0] = _eta[1] = _eta[2] = _eta[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixCountdown::end_of_step()
{
  const double wtime = MPI_Wtime();
  const double delta = wtime - _last_wall;
  const double tstep = static_cast<double>(update->laststep-update->firststep); 

  _last_wall = wtime;
  _times[_idx++] = delta;
  if (_ntimes < maxtimes) ++_ntimes;
  if (_idx == maxtimes) _idx = 0; 
  _complete = static_cast<double>(update->ntimestep-update->firststep)/tstep; 

  double tsum = 0.0;
  double wsum = 0.0;
  double w = 1.0;
  int i,j;

  // compute weighted average of time per timestep
  for (i=_idx; i < _idx+_ntimes; ++i) {
    j = (i < _ntimes) ? i : i - _ntimes;
    tsum += w*_times[j];
    wsum += w;
    w *= decay;
  }
  w = tsum/wsum * tstep * (1.0 - _complete) / static_cast<double>(nevery);
  int day = floor(w/86400.0);
  _eta[0] = static_cast<double>(day);
  w -= _eta[0]*86400.0;
  int hrs = floor(w/3600.0);
  _eta[1] = static_cast<double>(hrs);
  w -= _eta[1]*3600.0;
  int min = floor(w/60.0);
  _eta[2] = static_cast<double>(min);
  w -= _eta[2]*60.0;
  _eta[3] = w;

  if (_outfreq > 0) {
    if ((update->ntimestep % _outfreq == 0) && (comm->me == 0)) {
      const char fmt[] = " %6.2f%% complete. "
         "Estimated time to completion: %dd %02dh %02dm %04.2fs\n";
      if (screen) fprintf(screen,fmt,_complete*100.0,day,hrs,min,w);
      if (logfile) fprintf(logfile,fmt,_complete*100.0,day,hrs,min,w);
    }
  }
}

/* ----------------------------------------------------------------------
   return completion fraction
------------------------------------------------------------------------- */

double FixCountdown::compute_scalar()
{
  return _complete;
}

/* ----------------------------------------------------------------------
   return estimated time to completion component
------------------------------------------------------------------------- */

double FixCountdown::compute_vector(int id)
{
  return _eta[id];
}

