/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NPAIR_OMP_H
#define LMP_NPAIR_OMP_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "comm.h"
#include "fix_omp.h"
#include "modify.h"
#include "thr_data.h"
#include "timer.h"

namespace LAMMPS_NS {

// these macros hide some ugly and redundant OpenMP related stuff
#if defined(_OPENMP)

// get access to number of threads and per-thread data structures via FixOMP
#define NPAIR_OMP_INIT                 \
  const int nthreads = comm->nthreads; \
  omp_set_num_threads(nthreads); \
  const int ifix = modify->find_fix("package_omp")

// get thread id and then assign each thread a fixed chunk of atoms
#define NPAIR_OMP_SETUP(num)                                           \
  {                                                                    \
    const int tid = omp_get_thread_num();                              \
    const int idelta = 1 + num / nthreads;                             \
    const int ifrom = tid * idelta;                                    \
    const int ito = ((ifrom + idelta) > num) ? num : (ifrom + idelta); \
    FixOMP *fix = static_cast<FixOMP *>(modify->fix[ifix]);            \
    ThrData *thr = fix->get_thr(tid);                                  \
    thr->timer(Timer::START);

#define NPAIR_OMP_CLOSE     \
  thr->timer(Timer::NEIGH); \
  }

#else /* !defined(_OPENMP) */

#define NPAIR_OMP_INIT

#define NPAIR_OMP_SETUP(num) \
  const int tid = 0;         \
  const int ifrom = 0;       \
  const int ito = num

#define NPAIR_OMP_CLOSE

#endif

}    // namespace LAMMPS_NS

#endif
