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

#ifndef LMP_PAIR_OMP_H
#define LMP_PAIR_OMP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOMP : public Pair {

 protected:
  double *eng_vdwl_thr;         // per thread accumulated vdw energy
  double *eng_coul_thr;         // per thread accumulated coulomb energies
  double **virial_thr;          // per thread virial
  double **eatom_thr;		// per thread per atom energy
  double ***vatom_thr;		// per thread per atom virial

  int maxeatom_thr, maxvatom_thr;
  
 public:
  PairOMP(class LAMMPS *);
  virtual ~PairOMP();

  virtual double memory_usage();

 protected:
  // threading adapted versions of the ev_tally infrastructure.
  void ev_setup_thr(int, int);
  void ev_reduce_thr();
  void ev_tally_thr(int, int, int, int, double, double, double,
		    double, double, double, int);

  // reduce per thread forces into the first part of the force
  // array that is used for the non-threaded parts and reset
  // the temporary storage to 0.0.
  // this is in the header to be inlined.
  // need to post a barrier to wait until all threads are done
  // with computing forces. the reduction can be threaded as well.
  void force_reduce_thr(double **fall, int nall, int nthreads, int tid)
    {
#if defined(_OPENMP)
#pragma omp barrier
    {
      double **f;
      const int idelta = (nthreads > 1) ? 1 + nall/nthreads : nall;
      const int ifrom = tid*idelta;
      const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
      for (int n = 1; n < nthreads; ++n) {
	const int toffs = n*nall;
	f = fall + toffs;
	for (int m = ifrom; m < ito; ++m) {
	  fall[m][0] += f[m][0];
	  f[m][0] = 0.0;
	  fall[m][1] += f[m][1];
	  f[m][1] = 0.0;
	  fall[m][2] += f[m][2];
	  f[m][2] = 0.0;
	}
      }
      }
#endif
    };
};

}

#endif
