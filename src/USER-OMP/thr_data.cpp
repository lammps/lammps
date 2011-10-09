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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ThrData::ThrData(int tid) : _tid(tid)
{
  eatom = vatom = NULL;
  _maxeatom = _maxvatom = 0;
}

/* ---------------------------------------------------------------------- */

ThrData::~ThrData()
{
  if (eatom) delete[] eatom;
  if (vatom) delete[] vatom;
}

/* ---------------------------------------------------------------------- */

void ThrData::clear(int nall)
{
  eng_vdwl = eng_coul = eng_bond = 0.0;
  memset(virial,0,6*sizeof(double));
  if (eatom)
    memset(eatom,0,nall*sizeof(double));
  if (vatom)
    memset(vatom,0,6*nall*sizeof(double));
}

/* ---------------------------------------------------------------------- */

void ThrData::grow_arrays(int nmax, int eflag, int vflag)
{
  if (eflag) {
    if (_maxeatom < nmax) {
      _maxeatom = nmax;
      if (eatom) delete[] eatom;
      eatom = new double[_maxeatom];
    }
  }
  if (vflag) {
    if (_maxvatom < nmax) {
      _maxvatom = nmax;
      if (vatom) delete[] vatom;
      vatom = new double[6*_maxvatom];
    }
  }
}

/* ---------------------------------------------------------------------- */

#if 0
double ThrData::memory_usage_thr() 
{
  const int nthreads=lmp->comm->nthreads;

  double bytes = nthreads * (3 + 7) * sizeof(double);
  bytes += nthreads * maxeatom_thr * sizeof(double);
  bytes += nthreads * maxvatom_thr * 6 * sizeof(double);
  return bytes;
}
#endif
