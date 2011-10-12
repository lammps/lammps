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
   per-thread data management for LAMMPS
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "thr_data.h"

#include <string.h>
#include <stdio.h>

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

void ThrData::check_tid(int tid)
{
  if (tid != _tid)
    fprintf(stderr,"WARNING: external and internal tid mismatch %d != %d\n",tid,_tid);
}

/* ---------------------------------------------------------------------- */

void ThrData::clear()
{
  eng_vdwl=eng_coul=eng_bond=eng_angle=eng_dihed=eng_imprp=eng_kspce=0.0;
  memset(virial_pair,0,6*sizeof(double));
  memset(virial_bond,0,6*sizeof(double));
  memset(virial_angle,0,6*sizeof(double));
  memset(virial_dihed,0,6*sizeof(double));
  memset(virial_imprp,0,6*sizeof(double));
  memset(virial_kspce,0,6*sizeof(double));

  eatom = NULL;
  vatom = NULL;
}

/* ---------------------------------------------------------------------- */

double ThrData::memory_usage() 
{
  double bytes = (7 + 6*6) * sizeof(double);
  bytes += 2 * sizeof(double*);
  bytes += 4 * sizeof(int);

  return bytes;
}

/* additional helper functions */

// reduce per thread data into the first part of the data
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on
// multi-dimensional arrays like force stored in this order
// x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with writing to the array .
void LAMMPS_NS::data_reduce_thr(double *dall, int nall, int nthreads, int ndim, int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    const int nvals = ndim*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nvals) ? nvals : (ifrom + idelta);

    for (int m = ifrom; m < ito; ++m) {
      for (int n = 1; n < nthreads; ++n) {
	dall[m] += dall[n*nvals + m];
	dall[n*nvals + m] = 0.0;
      }
    }
  }
#else
  // NOOP in non-threaded execution.
  return;
#endif
}

