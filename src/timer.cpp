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

#include "timer.h"
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "comm.h"
#include "error.h"
#include "force.h"

#ifdef _WIN32
#include <windows.h>
#include <cstdint>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <ctime>

using namespace LAMMPS_NS;

// convert a timespec ([[HH:]MM:]SS) to seconds
// the strings "off" and "unlimited" result in -1;

static double timespec2seconds(char *timespec)
{
  double vals[3];
  char *num;
  int i = 0;

  // first handle allowed textual inputs
  if (strcmp(timespec,"off") == 0) return -1;
  if (strcmp(timespec,"unlimited") == 0) return -1;

  vals[0] = vals[1] = vals[2] = 0;

  num = strtok(timespec,":");
  while ((num != NULL) && (i < 3)) {
    vals[i] = atoi(num);
    ++i;
    num = strtok(NULL,":");
  }

  if (i == 3) return (vals[0]*60 + vals[1])*60 + vals[2];
  else if (i == 2) return vals[0]*60 + vals[1];
  else return vals[0];
}

// Return the CPU time for the current process in seconds very
// much in the same way as MPI_Wtime() returns the wall time.

static double CPU_Time()
{
  double rv = 0.0;

#ifdef _WIN32

  // from MSD docs.
  FILETIME ct,et,kt,ut;
  union { FILETIME ft; uint64_t ui; } cpu;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }

#else /* ! _WIN32 */

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }

#endif /* ! _WIN32 */

  return rv;
}

/* ---------------------------------------------------------------------- */

Timer::Timer(LAMMPS *lmp) : Pointers(lmp)
{
  _level = NORMAL;
  _sync  = OFF;
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
  double current_cpu=0.0, current_wall=0.0;

  if (_level > NORMAL) current_cpu = CPU_Time();
  current_wall = MPI_Wtime();

  if ((which > TOTAL) && (which < NUM_TIMER)) {
    const double delta_cpu = current_cpu - previous_cpu;
    const double delta_wall = current_wall - previous_wall;

    cpu_array[which]  += delta_cpu;
    wall_array[which] += delta_wall;
    cpu_array[ALL]    += delta_cpu;
    wall_array[ALL]   += delta_wall;
  }

  previous_cpu  = current_cpu;
  previous_wall = current_wall;

  if (which == RESET) {
    this->init();
    cpu_array[TOTAL] = current_cpu;
    wall_array[TOTAL] = current_wall;
  }

  if (_sync) {
    MPI_Barrier(world);
    if (_level > NORMAL) current_cpu = CPU_Time();
    current_wall = MPI_Wtime();

    cpu_array[SYNC]  += current_cpu - previous_cpu;
    wall_array[SYNC] += current_wall - previous_wall;
    previous_cpu  = current_cpu;
    previous_wall = current_wall;
  }
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start()
{
  double current_cpu=0.0, current_wall=0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = CPU_Time();
  current_wall = MPI_Wtime();

  cpu_array[TOTAL]  = current_cpu;
  wall_array[TOTAL] = current_wall;
  previous_cpu  = current_cpu;
  previous_wall = current_wall;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop()
{
  double current_cpu=0.0, current_wall=0.0;

  MPI_Barrier(world);

  if (_level < LOOP) return;

  current_cpu = CPU_Time();
  current_wall = MPI_Wtime();

  cpu_array[TOTAL]  = current_cpu - cpu_array[TOTAL];
  wall_array[TOTAL] = current_wall - wall_array[TOTAL];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu(enum ttype which)
{
  double current_cpu = CPU_Time();
  return (current_cpu - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(enum ttype which)
{
  if (_level == OFF) return 0.0;
  double current_wall = MPI_Wtime();
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
    const double d = MPI_Wtime() - timeout_start;
    // remaining timeout in seconds
    int s = _timeout - d;
    // remaining 1/100ths of seconds
    const int hs = 100*((_timeout - d) -  s);
    // breaking s down into second/minutes/hours
    const int seconds = s % 60;
    s  = (s - seconds) / 60;
    const int minutes = s % 60;
    const int hours = (s - minutes) / 60;
    fprintf(fp,"  Walltime left : %d:%02d:%02d.%02d\n",
            hours,minutes,seconds,hs);
  }
}

/* ---------------------------------------------------------------------- */

bool Timer::_check_timeout()
{
  double walltime = MPI_Wtime() - timeout_start;
  // broadcast time to insure all ranks act the same.
  MPI_Bcast(&walltime,1,MPI_DOUBLE,0,world);

  if (walltime < _timeout) {
    _nextcheck += _checkfreq;
    return false;
  } else {
    if (comm->me == 0)
      error->warning(FLERR,"Wall time limit reached");
    _timeout = 0.0;
    return true;
  }
}

/* ---------------------------------------------------------------------- */
double Timer::get_timeout_remain()
{
  return (_timeout < 0.0) ? 0.0 : _timeout + timeout_start - MPI_Wtime();
}

/* ----------------------------------------------------------------------
   modify parameters of the Timer class
------------------------------------------------------------------------- */
static const char *timer_style[] = { "off", "loop", "normal", "full" };
static const char *timer_mode[]  = { "nosync", "(dummy)", "sync" };
static const char  timer_fmt[]   = "New timer settings: style=%s  mode=%s  timeout=%s\n";

void Timer::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],timer_style[OFF])           == 0) {
      _level = OFF;
    } else if (strcmp(arg[iarg],timer_style[LOOP]) == 0) {
      _level = LOOP;
    } else if (strcmp(arg[iarg],timer_style[NORMAL]) == 0) {
      _level = NORMAL;
    } else if (strcmp(arg[iarg],timer_style[FULL])   == 0) {
      _level = FULL;
    } else if (strcmp(arg[iarg],timer_mode[OFF])     == 0) {
      _sync  = OFF;
    } else if (strcmp(arg[iarg],timer_mode[NORMAL])  == 0) {
      _sync  = NORMAL;
    } else if (strcmp(arg[iarg],"timeout") == 0) {
      ++iarg;
      if (iarg < narg) {
        _timeout = timespec2seconds(arg[iarg]);
      } else error->all(FLERR,"Illegal timers command");
    } else if (strcmp(arg[iarg],"every") == 0) {
      ++iarg;
      if (iarg < narg) {
        _checkfreq = force->inumeric(FLERR,arg[iarg]);
        if (_checkfreq <= 0)
          error->all(FLERR,"Illegal timers command");
      } else error->all(FLERR,"Illegal timers command");
    } else error->all(FLERR,"Illegal timers command");
    ++iarg;
  }

  timeout_start = MPI_Wtime();
  if (comm->me == 0) {

    // format timeout setting
    char timebuf[32];
    if (_timeout < 0) strcpy(timebuf,"off");
    else {
      time_t tv = _timeout;
      struct tm *tm = gmtime(&tv);
      strftime(timebuf,32,"%H:%M:%S",tm);
    }

    if (screen)
      fprintf(screen,timer_fmt,timer_style[_level],timer_mode[_sync],timebuf);
    if (logfile)
      fprintf(logfile,timer_fmt,timer_style[_level],timer_mode[_sync],timebuf);
  }
}
