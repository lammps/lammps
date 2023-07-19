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

#include "timer.h"

#include "comm.h"
#include "error.h"
#include "fmt/chrono.h"

#include <cstring>
#include <ctime>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(LAMMPS *_lmp) : Pointers(_lmp)
{
  _level = NORMAL;
  _sync = OFF;
  _timeout = -1;
  _s_timeout = -1;
  _checkfreq = 10;
  _nextcheck = -1;
  this->_stamp(RESET);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < NUM_TIMER; i++) {
    cpu_array[i] = 0.0;
    wall_array[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::_stamp(enum ttype which)
{
  double current_cpu = 0.0, current_wall = 0.0;

  if (_level > NORMAL) current_cpu = platform::cputime();
  current_wall = platform::walltime();

  if ((which > TOTAL) && (which < NUM_TIMER)) {
    const double delta_cpu = current_cpu - previous_cpu;
    const double delta_wall = current_wall - previous_wall;

    cpu_array[which] += delta_cpu;
    wall_array[which] += delta_wall;
    cpu_array[ALL] += delta_cpu;
    wall_array[ALL] += delta_wall;
  }

  previous_cpu = current_cpu;
  previous_wall = current_wall;

  if (which == RESET) {
    this->init();
    cpu_array[TOTAL] = current_cpu;
    wall_array[TOTAL] = current_wall;
  }

  if (_sync) {
    MPI_Barrier(world);
    if (_level > NORMAL) current_cpu = platform::cputime();
    current_wall = platform::walltime();

    cpu_array[SYNC] += current_cpu - previous_cpu;
    wall_array[SYNC] += current_wall - previous_wall;
    previous_cpu = current_cpu;
    previous_wall = current_wall;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start()
{
  double current_cpu = 0.0, current_wall = 0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = platform::cputime();
  current_wall = platform::walltime();

  cpu_array[TOTAL] = current_cpu;
  wall_array[TOTAL] = current_wall;
  previous_cpu = current_cpu;
  previous_wall = current_wall;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop()
{
  double current_cpu = 0.0, current_wall = 0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = platform::cputime();
  current_wall = platform::walltime();

  cpu_array[TOTAL] = current_cpu - cpu_array[TOTAL];
  wall_array[TOTAL] = current_wall - wall_array[TOTAL];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu(enum ttype which)
{
  double current_cpu = platform::cputime();
  return (current_cpu - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(enum ttype which)
{
  if (_level == OFF) return 0.0;
  double current_wall = platform::walltime();
  return (current_wall - wall_array[which]);
}

/* ---------------------------------------------------------------------- */

void Timer::set_wall(enum ttype which, double newtime)
{
  wall_array[which] = newtime;
}

/* ---------------------------------------------------------------------- */

void Timer::init_timeout()
{
  _s_timeout = _timeout;
  if (_timeout < 0)
    _nextcheck = -1;
  else
    _nextcheck = _checkfreq;
}

/* ---------------------------------------------------------------------- */

void Timer::print_timeout(FILE *fp)
{
  if (!fp) return;

  // format timeout setting
  if (_timeout > 0) {
    // time since init_timeout()
    const double d = platform::walltime() - timeout_start;
    // remaining timeout in seconds
    int s = _timeout - d;
    // remaining 1/100ths of seconds
    const int hs = 100 * ((_timeout - d) - s);
    // breaking s down into second/minutes/hours
    const int seconds = s % 60;
    s = (s - seconds) / 60;
    const int minutes = s % 60;
    const int hours = (s - minutes) / 60;
    fprintf(fp, "  Walltime left : %d:%02d:%02d.%02d\n", hours, minutes, seconds, hs);
  }
}

/* ---------------------------------------------------------------------- */

bool Timer::_check_timeout()
{
  double walltime = platform::walltime() - timeout_start;
  // broadcast time to ensure all ranks act the same.
  MPI_Bcast(&walltime, 1, MPI_DOUBLE, 0, world);

  if (walltime < _timeout) {
    _nextcheck += _checkfreq;
    return false;
  } else {
    if (comm->me == 0) error->warning(FLERR, "Wall time limit reached");
    _timeout = 0.0;
    return true;
  }
}

/* ---------------------------------------------------------------------- */
double Timer::get_timeout_remain()
{
  return (_timeout < 0.0) ? 0.0 : _timeout + timeout_start - platform::walltime();
}

/* ----------------------------------------------------------------------
   modify parameters of the Timer class
------------------------------------------------------------------------- */
static const char *timer_style[] = {"off", "loop", "normal", "full"};
static const char *timer_mode[] = {"nosync", "(dummy)", "sync"};

void Timer::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], timer_style[OFF]) == 0) {
      _level = OFF;
    } else if (strcmp(arg[iarg], timer_style[LOOP]) == 0) {
      _level = LOOP;
    } else if (strcmp(arg[iarg], timer_style[NORMAL]) == 0) {
      _level = NORMAL;
    } else if (strcmp(arg[iarg], timer_style[FULL]) == 0) {
      _level = FULL;
    } else if (strcmp(arg[iarg], timer_mode[OFF]) == 0) {
      _sync = OFF;
    } else if (strcmp(arg[iarg], timer_mode[NORMAL]) == 0) {
      _sync = NORMAL;
    } else if (strcmp(arg[iarg], "timeout") == 0) {
      ++iarg;
      if (iarg < narg) {
        _timeout = utils::timespec2seconds(arg[iarg]);
      } else
        error->all(FLERR, "Illegal timer command");
    } else if (strcmp(arg[iarg], "every") == 0) {
      ++iarg;
      if (iarg < narg) {
        _checkfreq = utils::inumeric(FLERR, arg[iarg], false, lmp);
        if (_checkfreq <= 0) error->all(FLERR, "Illegal timer command");
      } else
        error->all(FLERR, "Illegal timer command");
    } else
      error->all(FLERR, "Illegal timer command");
    ++iarg;
  }

  timeout_start = platform::walltime();
  if (comm->me == 0) {

    // format timeout setting
    std::string timeout = "off";
    if (_timeout >= 0) {
      std::tm tv = fmt::gmtime((std::time_t) _timeout);
      timeout = fmt::format("{:02d}:{:%M:%S}", tv.tm_yday * 24 + tv.tm_hour, tv);
    }

    utils::logmesg(lmp, "New timer settings: style={}  mode={}  timeout={}\n", timer_style[_level],
                   timer_mode[_sync], timeout);
  }
}
