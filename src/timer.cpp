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

#include "mpi.h"
#include "timer.h"
#include "memory.h"

// #define LMP_CLOCK_GETTIME 1

#ifdef LMP_CLOCK_GETTIME
#include <time.h>
#endif
#ifndef _WIN32
#include <sys/time.h>
#include <sys/resource.h>
#endif

using namespace LAMMPS_NS;

static double get_cpu_time()
{
#ifdef LMP_CLOCK_GETTIME
  struct timespec tp;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tp);
  double rv = (double) tp.tv_sec;
  rv += (double) tp.tv_nsec * 0.0000000001;
  return rv;
#else
#ifndef _WIN32
  struct rusage ru;
  double rv = 0.0;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }
  return rv;
#else  /* _WIN32 */
  // from MSD docs.
  FILETIME ct,et,kt,ut;
  double rv = 0.0;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    rv = (double) dwLowDateTime * 0.0000001;
    rv += (double) dwHighDateTime *429.4967296;
  }
  return rv;
);
#endif /* _WIN32 */
#endif
}

/* ---------------------------------------------------------------------- */
static double get_wall_time()
{
#ifdef LMP_CLOCK_GETTIME
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME,&tp);
  double rv = (double) tp.tv_sec;
  rv += (double) tp.tv_nsec * 0.0000000001;
  return rv;
#else
  return MPI_Wtime();
#endif
}

/* ---------------------------------------------------------------------- */

Timer::Timer(LAMMPS *lmp) : Pointers(lmp)
{
  memory->create(cpu_array,TIME_N,"timer:cpu_array");
  memory->create(wall_array,TIME_N,"timer:wall_array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->destroy(cpu_array);
  memory->destroy(wall_array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) cpu_array[i] = 0.0;
  for (int i = 0; i < TIME_N; i++) wall_array[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_cpu  = get_cpu_time();
  previous_wall = get_wall_time();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(enum ttype which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_cpu  = get_cpu_time();
  double current_wall = get_wall_time();
  cpu_array[which] += current_cpu - previous_cpu;
  wall_array[which] += current_wall - previous_wall;
  previous_cpu = current_cpu;
  previous_wall = current_wall;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(enum ttype which)
{
  MPI_Barrier(world);
  cpu_array[which] = get_cpu_time();
  wall_array[which] = get_wall_time();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(enum ttype which)
{
  MPI_Barrier(world);
  double current_cpu = get_cpu_time();
  double current_wall = get_wall_time();
  cpu_array[which] = current_cpu - cpu_array[which];
  wall_array[which] = current_wall - wall_array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu(enum ttype which)
{
  double current_cpu = get_cpu_time();
  return (current_cpu - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(enum ttype which)
{
  double current_wall = get_wall_time();
  return (current_wall - wall_array[which]);
}

/* ---------------------------------------------------------------------- */

void Timer::set_wall(enum ttype which, double newtime)
{
  wall_array[which] = newtime;
}

