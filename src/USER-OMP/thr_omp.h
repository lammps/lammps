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
#include "fix_omp.h"
#include "thr_data.h"

namespace LAMMPS_NS {

// forward declarations
class Pair;
class Bond;
class Angle;
class Dihedral;
class Improper;
class Kspace;

class ThrOMP {
 protected:
  LAMMPS *lmp; // reference to base lammps object.
  FixOMP *fix; // pointer to fix_omp;

  const int thr_style;

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

 protected:
  // extra ev_tally setup work for threaded styles
  void ev_setup_thr(int, int, int, double *, double **, ThrData *);

  ////////////////////////////////////////////////////////////////////////
  //  helper functions operating on data replicated for thread support  //
  ////////////////////////////////////////////////////////////////////////
  // compute global per thread virial contribution from per-thread force
  void virial_fdotr_compute_thr(double * const, const double * const * const, 
				const double * const * const,
				const int, const int, const int);

 private:
  // internal method to be used by multiple ev_setup_thr() methods
  void ev_setup_acc_thr(int, int, int, int, int, int);

 protected:
  // threading adapted versions of the ev_tally infrastructure
  // style specific versions (need access to style class flags)
  void ev_tally_thr(Pair *, int, int, int, int, double, double,
		    double, double, double, double, ThrData *);
  void ev_tally_xyz_thr(Pair *, int, int, int, int, double, double,
			double, double, double, double, double, double, int);
  void ev_tally3_thr(Pair *, int, int, int, double, double,
		     double *, double *, double *, double *, int);
  void ev_tally4_thr(Pair *, int, int, int, int, double, 
		     double *, double *, double *,
		     double *, double *, double *, int);
  void ev_tally_list_thr(Pair *, int, int *, double , double *, int);

  void ev_tally_thr(Dihedral *, int, int, int, int, int, int, double,
		    double *, double *, double *, double, double, double,
		    double, double, double, double, double, double, int);

  // style independent versions
  void v_tally2_thr(int, int, double, double *, int);
  void v_tally3_thr(int, int, int, double *, double *, double *, double *, int);
  void v_tally4_thr(int, int, int, int, double *, double *, double *,
		    double *, double *, double *, int);

};

}
#endif
