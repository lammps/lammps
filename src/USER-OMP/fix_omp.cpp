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
   OpenMP based threading support for LAMMPS
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "integrate.h"
#include "min.h"

#include "fix_omp.h"
#include "thr_data.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

using namespace LAMMPS_NS;

static int get_tid()
{
  int tid = 0;
#if defined(_OPENMP)
  tid = omp_get_thread_num();
#endif
  return tid;
}

/* ---------------------------------------------------------------------- */

FixOMP::FixOMP(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg),
  thr(NULL)
{
  if ((narg < 4) || (narg > 6)) error->all(FLERR,"Illegal fix OMP command");
  if (strcmp(arg[1],"all") != 0) error->all(FLERR,"Illegal fix OMP command");

  int nthreads = 1;
#if defined(_OPENMP)
  if (strcmp(arg[3],"*") == 0)
#pragma omp parallel default(none) shared(nthreads)
    nthreads = omp_get_num_threads();
  else
    nthreads = atoi(arg[3]);
#endif

  if (nthreads < 1)
    error->all(FLERR,"Illegal number of threads requested.");
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#endif
  comm->nthreads = nthreads;
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"  reset %d OpenMP thread(s) per MPI task\n", nthreads);
    if (logfile)
      fprintf(logfile,"  reset %d OpenMP thread(s) per MPI task\n", nthreads);
  }

#if 0 /* to be enabled when we can switch between half and full neighbor lists */
  _newton = true;
  if (narg > 5) {
    if (strcmp(arg[4],"newton") == 0) {
      if (strcmp(arg[5],"off") == 0)
	_newton = false;
      else if (strcmp(arg[5],"on") == 0)
	_newton = true;
      else
	error->all(FLERR,"Illegal fix OMP command");
    } else error->all(FLERR,"Illegal fix OMP command");
  }
#endif

  // allocate per thread accumulator manager class
  thr = new ThrData *[nthreads];
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    const int tid = get_tid();
    thr[tid] = new ThrData(tid);
  }
}

/* ---------------------------------------------------------------------- */

FixOMP::~FixOMP()
{
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    const int tid = get_tid();
    delete thr[tid];
  }
  delete[] thr;
}

/* ---------------------------------------------------------------------- */

int FixOMP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
#if 0
  mask |= PRE_FORCE_RESPA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
#endif
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOMP::init()
{
  fprintf(stderr,"%s:\n",__FUNCTION__);
  last_omp_style = THR_NONE;

#define CheckStyleForOMP(name,def)			\
  if (force->name) {					\
    int len = strlen(force->name ## _style);		\
    char *suffix = force->name ## _style + len - 4;	\
    if (strcmp(suffix,"/omp") == 0)			\
      last_omp_style = def;				\
  }

  // determine which is the last force style with OpenMP
  // support as this is the one that has to reduce the forces
  CheckStyleForOMP(pair,THR_PAIR);
  CheckStyleForOMP(bond,THR_BOND);
  CheckStyleForOMP(angle,THR_ANGLE);
  CheckStyleForOMP(dihedral,THR_DIHEDRAL);
  CheckStyleForOMP(improper,THR_IMPROPER);
  CheckStyleForOMP(kspace,THR_KSPACE);

#undef CheckStyleForOMP
}

/* ---------------------------------------------------------------------- */

void FixOMP::setup_pre_force(int vflag)
{
  fprintf(stderr,"%s: vflag=%d\n",__FUNCTION__, vflag);
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixOMP::setup(int vflag)
{
  fprintf(stderr,"%s: vflag=%d\n",__FUNCTION__, vflag);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

// adjust size and clear out per thread accumulator arrays
void FixOMP::pre_force(int)
{
  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  
  double * const f = &(atom->f[0][0]);
  
  memset(f,0,3*nthreads*nall*sizeof(double));
  
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    const int tid = get_tid();
    thr[tid]->check_tid(tid);
    thr[tid]->clear();
  }
}

/* ---------------------------------------------------------------------- */

void FixOMP::post_force(int vflag)
{
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  double * const f = &(atom->f[0][0]);

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    const int tid = get_tid();
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
    if (atom->torque)
      data_reduce_thr(&(atom->torque[0][0]), nall, nthreads, 3, tid);
  }
}


/* ---------------------------------------------------------------------- */

double FixOMP::memory_usage()
{
  double bytes = comm->nthreads * sizeof(ThrData *);
  bytes += comm->nthreads * thr[0]->memory_usage();
  
  return bytes;
}
