/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_THR_OMP_H
#define LMP_THR_OMP_H

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "error.h"
#include "fix_omp.h"    // IWYU pragma: export
#include "pointers.h"
#include "thr_data.h"    // IWYU pragma: export

#include <atomic>
#include <string>

namespace LAMMPS_NS {

// forward declarations
class Pair;
class Bond;
class Angle;
class Dihedral;
class Improper;

class ThrOMP {

 protected:
  LAMMPS *lmp;    // reference to base lammps object.
  FixOMP *fix;    // pointer to fix_omp;

  int thr_style;
  // use a C++11 atomic, since OpenMP 3.1 atomic reads are not available with all compilers
  std::atomic<int> thr_error;
  int thr_errline;
  const char *thr_errfile;
  std::string thr_errmsg;

 public:
  ThrOMP(LAMMPS *, int);
  // clang-format off
  // Cannot use = default here due to broken GCC 8 on RHEL 8
  virtual ~ThrOMP() noexcept(false) {} // NOLINT
  // clang-format on

  double memory_usage_thr();

  inline void sync_threads(){
#if defined(_OPENMP)
#pragma omp barrier
#endif
      {;
}
};    // namespace LAMMPS_NS

enum {
  THR_NONE = 0,
  THR_PAIR = 1,
  THR_BOND = 1 << 1,
  THR_ANGLE = 1 << 2,
  THR_DIHEDRAL = 1 << 3,
  THR_IMPROPER = 1 << 4,
  THR_KSPACE = 1 << 5,
  THR_CHARMM = 1 << 6, /*THR_PROXY=1<<7,THR_HYBRID=1<<8, */
  THR_FIX = 1 << 9,
  THR_INTGR = 1 << 10
};

protected:
// extra ev_tally setup work for threaded styles
void ev_setup_thr(int, int, int, double *, double **, double **, ThrData *);

// compute global per thread virial contribution from per-thread force
void virial_fdotr_compute_thr(double *const, const double *const *const, const double *const *const,
                              const int, const int, const int);

// reduce per thread data as needed
void reduce_thr(void *const style, const int eflag, const int vflag, ThrData *const thr);

// thread safe variant error abort support.
// signals an error condition in any thread by making
// thr_error > 0, if condition "cond" is true, and records
// the location and a copy of the message of the first error.
// returns true if an error was signaled by any thread,
// otherwise false. use return value to jump/return to the
// end of the threaded region and call error_thr() after
// the threaded region to stop with the recorded error.

bool check_error_thr(const bool cond, const int tid, const char *fname, const int line,
                     const std::string &errmsg)
{
  return check_error_thr(cond, tid, fname, line, errmsg.c_str());
}

bool check_error_thr(const bool cond, const int /*tid*/, const char *fname, const int line,
                     const char *errmsg)
{
  if (cond) {
    // only the first thread to signal an error records its location and message.
    // they are read by error_thr() after the threaded region has ended.
    if (thr_error.fetch_add(1, std::memory_order_relaxed) == 0) {
      thr_errfile = fname;
      thr_errline = line;
      thr_errmsg = errmsg;
    }
    return true;
  }
  return thr_error.load(std::memory_order_relaxed) > 0;
}

// stop with the first error that was signaled by check_error_thr().
// must be called outside of threaded regions.

void error_thr()
{
  if (thr_error.load() > 0) {
    thr_error = 0;
    lmp->error->one(thr_errfile, thr_errline, thr_errmsg);
  }
}

protected:
// threading adapted versions of the ev_tally infrastructure
// style specific versions (need access to style class flags)

// Pair
void e_tally_thr(Pair *const, const int, const int, const int, const int, const double,
                 const double, ThrData *const);
void v_tally_thr(Pair *const, const int, const int, const int, const int, const double *const,
                 ThrData *const);

void ev_tally_thr(Pair *const, const int, const int, const int, const int, const double,
                  const double, const double, const double, const double, const double,
                  ThrData *const);
void ev_tally_full_thr(Pair *const, const int, const double, const double, const double,
                       const double, const double, const double, ThrData *const);
void ev_tally_xyz_thr(Pair *const, const int, const int, const int, const int, const double,
                      const double, const double, const double, const double, const double,
                      const double, const double, ThrData *const);
void ev_tally_xyz_full_thr(Pair *const, const int, const double, const double, const double,
                           const double, const double, const double, const double, const double,
                           ThrData *const);
void v_tally2_thr(Pair *const, const int, const int, const double, const double *const,
                  ThrData *const);
void v_tally2_newton_thr(Pair *const, const int, const double *const, const double *const,
                         ThrData *const);
void ev_tally3_thr(Pair *const, const int, const int, const int, const double, const double,
                   const double *const, const double *const, const double *const,
                   const double *const, ThrData *const);
void v_tally3_thr(Pair *const, const int, const int, const int, const double *const,
                  const double *const, const double *const, const double *const, ThrData *const);
void ev_tally4_thr(Pair *const, const int, const int, const int, const int, const double,
                   const double *const, const double *const, const double *const,
                   const double *const, const double *const, const double *const, ThrData *const);
void v_tally4_thr(Pair *const, const int, const int, const int, const int, const double *const,
                  const double *const, const double *const, const double *const,
                  const double *const, const double *const, ThrData *const);

// Bond
void ev_tally_thr(Bond *const, const int, const int, const int, const int, const double,
                  const double, const double, const double, const double, ThrData *const);

// Angle
void ev_tally_thr(Angle *const, const int, const int, const int, const int, const int, const double,
                  const double *const, const double *const, const double, const double,
                  const double, const double, const double, const double, ThrData *const thr);
void ev_tally13_thr(Angle *const, const int, const int, const int, const int, const double,
                    const double, const double, const double, const double, ThrData *const thr);

// Dihedral
void ev_tally_thr(Dihedral *const, const int, const int, const int, const int, const int, const int,
                  const double, const double *const, const double *const, const double *const,
                  const double, const double, const double, const double, const double,
                  const double, const double, const double, const double, ThrData *const);

// Improper
void ev_tally_thr(Improper *const, const int, const int, const int, const int, const int, const int,
                  const double, const double *const, const double *const, const double *const,
                  const double, const double, const double, const double, const double,
                  const double, const double, const double, const double, ThrData *const);

// style independent versions
void ev_tally_list_thr(Pair *const, const int, const int *const, const double *const, const double,
                       const double, ThrData *const);
};

#if defined(_OPENMP)
#define OMP_ONLY_PARAM(x) x
#else
#define OMP_ONLY_PARAM(x)
#endif

// set loop range thread id, and force array offset for threaded runs.
inline void loop_setup_thr(int &ifrom, int &ito, int &tid, int inum, int OMP_ONLY_PARAM(nthreads))
{
#if defined(_OPENMP)
  tid = omp_get_thread_num();

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + inum / nthreads;
  ifrom = tid * idelta;
  ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
  tid = 0;
  ifrom = 0;
  ito = inum;
#endif
}
#undef OMP_ONLY_PARAM

// helpful definitions to help compilers optimizing code better

using dbl3_t = struct _dbl3_t {
  double x, y, z;
};
using dbl4_t = struct _dbl4_t {
  double x, y, z, w;
};
using int3_t = struct _int3_t {
  int a, b, t;
};
using int4_t = struct _int4_t {
  int a, b, c, t;
};
using int5_t = struct _int5_t {
  int a, b, c, d, t;
};
}    // namespace LAMMPS_NS

#endif
