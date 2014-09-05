/* -*- c++ -*- -------------------------------------------------------------
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

#ifndef LMP_THR_OMP_H
#define LMP_THR_OMP_H

#include "pointers.h"
#include "error.h"
#include "fix_omp.h"
#include "thr_data.h"

namespace LAMMPS_NS {

// forward declarations
class Pair;
class Bond;
class Angle;
class Dihedral;
class Improper;
class KSpace;
class Fix;

class ThrOMP {

 protected:
  LAMMPS *lmp; // reference to base lammps object.
  FixOMP *fix; // pointer to fix_omp;

  const int thr_style;
  int thr_error;

 public:
  ThrOMP(LAMMPS *, int);
  virtual ~ThrOMP();

  double memory_usage_thr();

  inline void sync_threads() {
#if defined(_OPENMP)
#pragma omp barrier
#endif
      { ; }
    };

  enum {THR_NONE=0,THR_PAIR=1,THR_BOND=1<<1,THR_ANGLE=1<<2,
        THR_DIHEDRAL=1<<3,THR_IMPROPER=1<<4,THR_KSPACE=1<<5,
        THR_CHARMM=1<<6, /*THR_PROXY=1<<7,THR_HYBRID=1<<8, */
        THR_FIX=1<<9,THR_INTGR=1<<10};

 protected:
  // extra ev_tally setup work for threaded styles
  void ev_setup_thr(int, int, int, double *, double **, ThrData *);

  // compute global per thread virial contribution from per-thread force
  void virial_fdotr_compute_thr(double * const, const double * const * const,
                                const double * const * const,
                                const int, const int, const int);

  // reduce per thread data as needed
  void reduce_thr(void * const style, const int eflag, const int vflag,
                  ThrData * const thr);

  // thread safe variant error abort support.
  // signals an error condition in any thread by making
  // thr_error > 0, if condition "cond" is true.
  // will abort from thread 0 if thr_error is > 0
  // otherwise return true.
  // returns false if no error found on any thread.
  // use return value to jump/return to end of threaded region.

  bool check_error_thr(const bool cond, const int tid, const char *fname,
                       const int line, const char *errmsg) {
    if (cond) {
#if defined(_OPENMP)
#pragma omp atomic
      ++thr_error;
#endif
      if (tid > 0) return true;
      else lmp->error->one(fname,line,errmsg);
    } else {
      if (thr_error > 0) {
        if (tid == 0) lmp->error->one(fname,line,errmsg);
        else return true;
      } else return false;
    }
    return false;
  };

 protected:

  // threading adapted versions of the ev_tally infrastructure
  // style specific versions (need access to style class flags)

  // Pair
  void e_tally_thr(Pair * const, const int, const int, const int,
                   const int, const double, const double, ThrData * const);
  void v_tally_thr(Pair * const, const int, const int, const int,
                   const int, const double * const, ThrData * const);

  void ev_tally_thr(Pair * const, const int, const int, const int, const int,
                    const double, const double, const double, const double,
                    const double, const double, ThrData * const);
  void ev_tally_xyz_thr(Pair * const, const int, const int, const int,
                        const int, const double, const double, const double,
                        const double, const double, const double,
                        const double, const double, ThrData * const);
  void ev_tally3_thr(Pair * const, const int, const int, const int, const double,
                     const double, const double * const, const double * const,
                     const double * const, const double * const, ThrData * const);
  void ev_tally4_thr(Pair * const, const int, const int, const int, const int,
                     const double, const double * const, const double * const,
                     const double * const, const double * const, const double * const,
                     const double * const, ThrData * const);

  // Bond
  void ev_tally_thr(Bond * const, const int, const int, const int, const int,
                    const double, const double, const double, const double,
                    const double, ThrData * const);

  // Angle
  void ev_tally_thr(Angle * const, const int, const int, const int, const int, const int,
                    const double, const double * const, const double * const,
                    const double, const double, const double, const double, const double,
                    const double, ThrData * const thr);
  void ev_tally13_thr(Angle * const, const int, const int, const int, const int,
                      const double, const double, const double, const double,
                      const double, ThrData * const thr);

  // Dihedral
  void ev_tally_thr(Dihedral * const, const int, const int, const int, const int, const int,
                    const int, const double, const double * const, const double * const,
                    const double * const, const double, const double, const double,
                    const double, const double, const double, const double, const double,
                    const double, ThrData * const);

  // Improper
  void ev_tally_thr(Improper * const, const int, const int, const int, const int, const int,
                    const int, const double, const double * const, const double * const,
                    const double * const, const double, const double, const double,
                    const double, const double, const double, const double, const double,
                    const double, ThrData * const);

  // style independent versions
  void v_tally2_thr(const int, const int, const double, const double * const, ThrData * const);
  void v_tally3_thr(const int, const int, const int, const double * const, const double * const,
                    const double * const, const double * const, ThrData * const);
  void v_tally4_thr(const int, const int, const int, const int, const double * const,
                    const double * const, const double * const, const double * const,
                    const double * const, const double * const, ThrData * const);
  void ev_tally_list_thr(Pair * const, const int, const int * const,
                         const double * const, const double, const double,
                         ThrData * const);

};

// set loop range thread id, and force array offset for threaded runs.
static inline void loop_setup_thr(int &ifrom, int &ito, int &tid,
                                  int inum, int nthreads)
{
#if defined(_OPENMP)
  tid = omp_get_thread_num();

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + inum/nthreads;
  ifrom = tid*idelta;
  ito   = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
  tid = 0;
  ifrom = 0;
  ito = inum;
#endif
}

// helpful definitions to help compilers optimizing code better

typedef struct { double x,y,z;   } dbl3_t;
typedef struct { double x,y,z,w; } dbl4_t;
typedef struct { int a,b,t;      } int3_t;
typedef struct { int a,b,c,t;    } int4_t;
typedef struct { int a,b,c,d,t;  } int5_t;

}

#endif
