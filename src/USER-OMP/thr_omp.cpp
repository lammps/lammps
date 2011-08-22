/* -------------------------------------------------------------------------
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
   OpenMP based threading support for LAMMPS
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "thr_omp.h"

#include "memory.h"
#include "comm.h"

#include <stdio.h>

#if defined(_OPENMP)      
#include <omp.h>
#endif

using namespace LAMMPS_NS;

ThrOMP::ThrOMP(LAMMPS *lmp, int style) : Pointers (lmp), thr_style(style)
{
  maxeatom_thr = maxvatom_thr = 0;
  eng_vdwl_thr = eng_coul_thr = eng_bond_thr = NULL;
  virial_thr = eatom_thr = NULL;
  vatom_thr = NULL;
}

ThrOMP::~ThrOMP() 
{
  memory->destroy(eng_vdwl_thr);
  memory->destroy(eng_coul_thr);
  memory->destroy(eng_bond_thr);
  memory->destroy(virial_thr);
  memory->destroy(eatom_thr);
  memory->destroy(vatom_thr);
}

// set loop range for, thread id, and force array offset for threaded runs.
double **ThrOMP::loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
				const int inum, const int nall, const int nthreads)
{
#if defined(_OPENMP)
  if (nthreads > 1) {
    tid = omp_get_thread_num();

    // each thread works on a fixed chunk of atoms.
    const int idelta = 1 + inum/nthreads;
    ifrom = tid*idelta;
    ito   = ifrom + idelta;
    if (ito > inum)
      ito = inum;

    return f + nall*tid;

  } else {
#endif
    tid = 0;
    ifrom = 0;
    ito = inum;
    return f;
#if defined(_OPENMP)
  }
#endif
}

#if 0
  // threading adapted versions of the ev_tally infrastructure.
  void ev_setup_thr(int, int);
  void ev_reduce_thr();

  void ev_tally_thr(int, int, int, int, double, double, double,
		    double, double, double, int);
  void ev_tally_xyz_thr(int, int, int, int, double, double,
			double, double, double, double, double, double, int);
  void ev_tally3_thr(int, int, int, double, double,
		     double *, double *, double *, double *, int);
  void ev_tally4_thr(int, int, int, int, double,
		     double *, double *, double *, double *, double *, double *, int);
  void ev_tally_list_thr(int, int *, double, double *, int);
  void v_tally2_thr(int, int, double, double *, int);
  void v_tally3_thr(int, int, int, double *, double *, double *, double *, int);
  void v_tally4_thr(int, int, int, int, double *, double *, double *,
                    double *, double *, double *, int);
#endif


// reduce per thread forces into the first part of the force
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0.
// this is in the header to be inlined.
// need to post a barrier to wait until all threads are done
// with computing forces. the reduction can be threaded as well.
void ThrOMP::force_reduce_thr(double **fall, const int nall,
			      const int nthreads, const int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    double **f;
    const int idelta = 1 + nall/nthreads;
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
#else
  // NOOP in non-threaded execution.
  return;
#endif
}

#if 0
  // reduce per thread density into the first part of the rho
  // array that is used for the non-threaded parts. for use with EAM.
  // this is in the header to be inlined.
  // we need to post a barrier to wait until all threads are done.
  // the reduction can be threaded as well.
  void rho_reduce_thr(double *rho, const int nmax, const int nrange, 
		      const int nthreads, const int tid)
    {
#if defined(_OPENMP)
      // NOOP in non-threaded execution.
      if (nthreads == 1) return;
#pragma omp barrier
      {
	double *rho_thr;
	const int idelta = 1 + nrange/nthreads;
	const int ifrom = tid*idelta;
	const int ito   = ((ifrom + idelta) > nrange) ? nrange : (ifrom + idelta);
	for (int n = 1; n < nthreads; ++n) {
	  const int toffs = n*nmax;
	  rho_thr = rho + toffs;
	  for (int m = ifrom; m < ito; ++m)
	    rho[m] += rho_thr[m];
	}
      }
#else
      // NOOP in non-threaded execution.
      return;
#endif
    };
#endif

#if 0
// dihedral versions
  void ev_setup_thr(int, int);
  void ev_reduce_thr();
  void ev_tally_thr(int, int, int, int, int, int, double,
		double *, double *, double *, double, double, double,
		double, double, double, double, double, double, int);
#endif

