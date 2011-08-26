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

#ifndef _WIN32
#include "sys/time.h"
#include "sys/resource.h"
#endif

using namespace LAMMPS_NS;

static double get_cpu_time()
{
#ifndef _WIN32
  struct rusage ru;
  double rv = 0.0;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }
  return rv;
#endif
  return 0.0;
}

/* ---------------------------------------------------------------------- */
static double get_sys_time()
{
#ifndef _WIN32
  struct rusage ru;
  double rv = 0.0;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_stime.tv_sec;
    rv += (double) ru.ru_stime.tv_usec * 0.000001;
  }
  return rv;
#endif
  return 0.0;
}

/* ---------------------------------------------------------------------- */
static double get_wall_time()
{
  return MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

Timer::Timer(LAMMPS *lmp) : Pointers(lmp)
{
  memory->create(cpu_array,TIME_N,"timer:cpu_array");
  memory->create(sys_array,TIME_N,"timer:sys_array");
  memory->create(wall_array,TIME_N,"timer:wall_array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->destroy(cpu_array);
  memory->destroy(sys_array);
  memory->destroy(wall_array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) cpu_array[i] = 0.0;
  for (int i = 0; i < TIME_N; i++) sys_array[i] = 0.0;
  for (int i = 0; i < TIME_N; i++) wall_array[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_cpu  = get_cpu_time();
  previous_sys  = get_sys_time();
  previous_wall = get_wall_time();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(enum ttype which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_cpu  = get_cpu_time();
  double current_sys  = get_sys_time();
  double current_wall = get_wall_time();
  cpu_array[which] += current_cpu - previous_cpu;
  sys_array[which] += current_sys - previous_sys;
  wall_array[which] += current_wall - previous_wall;
  previous_cpu = current_cpu;
  previous_sys = current_sys;
  previous_wall = current_wall;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(enum ttype which)
{
  MPI_Barrier(world);
  cpu_array[which] = get_cpu_time();
  sys_array[which] = get_sys_time();
  wall_array[which] = get_wall_time();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(enum ttype which)
{
  MPI_Barrier(world);
  double current_cpu = get_cpu_time();
  double current_sys = get_sys_time();
  double current_wall = get_wall_time();
  cpu_array[which] = current_cpu - cpu_array[which];
  sys_array[which] = current_sys - sys_array[which];
  wall_array[which] = current_wall - wall_array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::cpu(enum ttype which)
{
  double current_cpu = get_cpu_time();
  return (current_cpu - cpu_array[which]);
}

/* ---------------------------------------------------------------------- */

double Timer::sys(enum ttype which)
{
  double current_sys = get_sys_time();
  return (current_sys - sys_array[which]);
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

