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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pppm_tip4p_proxy.h"

#include <stdio.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PPPMTIP4PProxy::PPPMTIP4PProxy(LAMMPS *lmp, int narg, char **arg) :
  PPPMTIP4P(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE|THR_PROXY) { need_setup=1; }

/* ---------------------------------------------------------------------- */

void PPPMTIP4PProxy::setup_proxy()
{
  if (need_setup)
    PPPMTIP4P::setup();

  need_setup=0;
}

/* ---------------------------------------------------------------------- */

void PPPMTIP4PProxy::compute(int eflag, int vflag)
{
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    ThrData *thr = fix->get_thr(tid);
    reduce_thr(this, eflag,vflag,thr);
  }
}

/* ---------------------------------------------------------------------- */

void PPPMTIP4PProxy::compute_proxy(int eflag, int vflag)
{
  setup_proxy();
  PPPMTIP4P::compute(eflag,vflag);
}

