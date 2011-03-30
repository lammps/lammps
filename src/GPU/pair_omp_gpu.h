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
   Modified by Mike for use with GPU library
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_OMP_GPU_H
#define LMP_PAIR_OMP_GPU_H

#include "pair.h"
#include <cassert>

#if defined(_OPENMP)      
#include "atom.h"
#include "memory.h"
#include <omp.h>
#endif

namespace LAMMPS_NS {

#if defined(_OPENMP)

class PairOMPGPU : protected Pointers {

 protected:
  double *eng_vdwl_thr;         // per thread accumulated vdw energy
  double *eng_coul_thr;         // per thread accumulated coulomb energies
  double **virial_thr;          // per thread virial
  double **eatom_thr;		// per thread per atom energy
  double ***vatom_thr;		// per thread per atom virial
  double **f_thr;               // per thread force (for thread id>0)

  int maxeatom_thr, maxvatom_thr;
  int _nthreads, _nmax;

  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;

 public:
  PairOMPGPU(LAMMPS *);
  ~PairOMPGPU();

  void mem_free();
  void init_style();
  double memory_usage();
 
  // set loop range for, thread id, and force array offset for threaded runs.
  double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
			  const int start, const int inum, const int nall)
    {
      if (_nthreads > 1) {
	tid = omp_get_thread_num();
	if (nall > _nmax) {
#pragma omp master
	  memory->grow(f_thr,atom->nmax*(_nthreads-1),3,"pair:f_thr");
#pragma omp barrier
	  if (tid == 0)
	    _nmax = atom->nmax;
	}

	// each thread works on a fixed chunk of atoms.
	const int idelta = 1 + (inum - start)/_nthreads;
	ifrom = tid * idelta + start;
	ito   = ifrom + idelta;
	if (ito > inum)
	  ito = inum;

	// zero per thread force array
	// keep thread memory access same as for force accumulation
	if (tid > 0) {
	  double **f_zero = f_thr + nall * (tid - 1);
	  for (int i = 0; i < nall; i++) {
	    f_zero[i][0] = 0.0;
	    f_zero[i][1] = 0.0;
	    f_zero[i][2] = 0.0;
	  }
	  return f_zero;
	}
	return f;
      } else {
	tid = 0;
	ifrom = start;
	ito = inum;
	return f;
      }
    };

  // set loop range for, thread id, and force array offset for threaded runs.
  double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
			  const int inum, const int nall)
  { return loop_setup_thr(f,ifrom,ito,tid,0,inum,nall); }

  // loop setup with full neighbor list and nonzero starting index
  double **loop_setup_thr_full(double **f, int &ifrom, int &ito, int &tid,
			       const int start, const int inum, const int nall)
    {
      if (_nthreads > 1) {
	tid = omp_get_thread_num();

	// each thread works on a fixed chunk of atoms.
	const int idelta = 1 + (inum - start) / _nthreads;
	ifrom = start + tid * idelta;
	ito   = ifrom + idelta;
	if (ito > inum)
	  ito = inum;

	return f;

      } else {
	tid = 0;
	ifrom = start;
	ito = inum;
	return f;
      }
    };

  // threading adapted versions of the ev_tally infrastructure.
  void ev_setup_thr(int, int, int, int, int, int, int, int);
  void ev_reduce_thr(Pair &);
  void ev_tally_thr(int, int, int, int, double, double, double,
		    double, double, double, int);
  void ev_tally_full_thr(int, double, double, double,
			 double, double, double, int);
  void ev_tally_xyz_thr(int, int, int, int, double, double,
			double, double, double, double, double, double, int);
  void ev_tally3_thr(int, int, int, double, double,
		     double *, double *, double *, double *, int, double);
  void ev_tally4_thr(int, int, int, int, double,
		     double *, double *, double *, double *, double *, 
		     double *, int);
  void ev_tally_list_thr(int, int *, double, double *, int);
  void v_tally2_thr(int, int, double, double *, int);
  void v_tally3_thr(int, int, int, double *, double *, double *, double *, int,
                    double);
  void v_tally4_thr(int, int, int, int, double *, double *, double *,
                    double *, double *, double *, int);

  // reduce per thread forces into the first part of the force
  // array that is used for the non-threaded parts and reset
  // the temporary storage to 0.0.
  // this is in the header to be inlined.
  // need to post a barrier to wait until all threads are done
  // with computing forces. the reduction can be threaded as well.
  void force_reduce_thr(double **fall, const int nall, const int tid)
    {
      // NOOP in non-threaded execution.
      if (_nthreads == 1) return;
#pragma omp barrier
      {
	double **f;
	const int idelta = 1 + nall/_nthreads;
	const int ifrom = tid*idelta;
	const int ito   = ((ifrom + idelta) > nall) ? nall : (ifrom + idelta);
	for (int n = 1; n < _nthreads; ++n) {
	  const int toffs = (n-1)*nall;
	  f = f_thr + toffs;
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
    };

  // reduce per thread density into the first part of the rho
  // array that is used for the non-threaded parts. for use with EAM.
  // this is in the header to be inlined.
  // we need to post a barrier to wait until all threads are done.
  // the reduction can be threaded as well.
  void rho_reduce_thr(double *rho, const int nmax, const int nrange, 
		      const int tid)
    {
      // NOOP in non-threaded execution.
      if (_nthreads == 1) return;
#pragma omp barrier
      {
	double *rho_thr;
	const int idelta = 1 + nrange/_nthreads;
	const int ifrom = tid*idelta;
	const int ito   = ((ifrom + idelta) > nrange) ? nrange : (ifrom + idelta);
	for (int n = 1; n < _nthreads; ++n) {
	  const int toffs = n*nmax;
	  rho_thr = rho + toffs;
	  for (int m = ifrom; m < ito; ++m)
	    rho[m] += rho_thr[m];
	}
      }
    };

};

#else
  
class PairOMPGPU {

 public:
  inline PairOMPGPU(LAMMPS *) {};
  virtual inline ~PairOMPGPU() {};

  inline void init_style() {}
  inline double memory_usage() { return 0.0; }

  inline void ev_setup_thr(int, int, int, int, int, int, int, int) {}
  inline void ev_reduce_thr(Pair &) {}
  inline double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
				 const int start, const int inum, 
				 const int nall)
    {
      ifrom = start;
      ito = inum;
      return f;
    };

  inline double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
				 const int inum, const int nall)
  { return loop_setup_thr(f,ifrom,ito,tid,0,inum,nall); }

  // loop setup with full neighbor list and nonzero starting index
  double **loop_setup_thr_full(double **f, int &ifrom, int &ito, int &tid,
			       const int start, const int inum, const int nall)
    {
      ifrom = start;
      ito = inum;
      return f;
    };

  inline void force_reduce_thr(double **fall, const int nall, const int tid) {};
  inline void rho_reduce_thr(double *rho, const int nmax, const int nrange, 
			     const int tid) {};

};

#endif

}

#endif

