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

#include "pppm_proxy.h"

#include <stdio.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PPPMProxy::PPPMProxy(LAMMPS *lmp, int narg, char **arg) :
  PPPM(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE|THR_PROXY) { need_setup=1; }

/* ---------------------------------------------------------------------- */

void PPPMProxy::setup_proxy()
{
  if (need_setup)
    PPPM::setup();

  need_setup = 0;
}

/* ---------------------------------------------------------------------- */

void PPPMProxy::compute(int eflag, int vflag)
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

void PPPMProxy::compute_proxy(int eflag, int vflag)
{
  setup_proxy();
  PPPM::compute(eflag,vflag);
}

