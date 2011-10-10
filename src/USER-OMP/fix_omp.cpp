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

#include "fix_omp.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "integrate.h"
#include "min.h"

#include "thr_data.h"

#include <string.h>
#include <stdio.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

static int get_tid()
{
  int tid = 0;
#if defined(_OPENMP)
  int tid = omp_get_thread_num();
#endif
  return tid;
}

/* ---------------------------------------------------------------------- */

FixOMP::FixOMP(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg),
  thr(NULL),_ndata(0), _nall(0), _nmax(0)
{
  if ((narg < 4) || (narg > 6)) error->all(FLERR,"Illegal fix OMP command");
  if (strcmp(arg[1],"all") != 0) error->all(FLERR,"Illegal fix OMP command");

  int nthreads = 1;
#if defined(_OPENMP)
  if (strcmp(arg[3],"*") == 0)
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
#pragma omp parallel default(none) shared(thr)
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
#pragma omp parallel default(none) shared(thr)
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
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= PRE_FORCE_RESPA;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;
#endif
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOMP::grow_arrays(int nmax)
{
  fprintf(stderr,"%s, nmax=%d->%d\n", __FUNCTION__, _nmax, nmax);
  _nmax=nmax;
}

/* ---------------------------------------------------------------------- */

void FixOMP::setup_pre_force(int vflag)
{
  fprintf(stderr,"%s, vflag=%d\n", __FUNCTION__, vflag);
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixOMP::setup_post_force(int vflag)
{
  fprintf(stderr,"%s, vflag=%d\n", __FUNCTION__, vflag);
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

// adjust size and clear out per thread accumulator arrays
void FixOMP::pre_force(int _vflag)
{
  fprintf(stderr,"%s, vflag=%d\n", __FUNCTION__, _vflag);
  _nall = atom->nlocal + atom->nghost;
  const int nall = _nall;
  const int nmax = _nmax;

  int accflags = ThrData::THR_NONE;
  if (force->pair)     accflags |= ThrData::THR_PAIR;
  if (force->bond)     accflags |= ThrData::THR_BOND;
  if (force->angle)    accflags |= ThrData::THR_ANGLE;
  if (force->dihedral) accflags |= ThrData::THR_DIHEDRAL;
  if (force->improper) accflags |= ThrData::THR_IMPROPER;
  if (force->kspace)   accflags |= ThrData::THR_KSPACE;

  int eflag,vflag;
  if (update->integrate) {
    eflag = update->integrate->eflag;
    vflag = update->integrate->vflag;
  } else if (update->minimize) {
    eflag = update->minimize->eflag;
    vflag = update->minimize->vflag;
  }
  if (eflag & 1) accflags |= ThrData::THR_ENERGY;
  if (eflag & 2) accflags |= ThrData::THR_EATOM;
  if (vflag & 3) accflags |= ThrData::THR_VIRIAL;
  if (vflag & 2) accflags |= ThrData::THR_VFDOTR;
  if (vflag & 4) accflags |= ThrData::THR_VATOM;
  
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(thr,accflags)
#endif
  {
    const int tid = get_tid();
    thr[tid]->set_accflags(accflags);
    thr[tid]->grow_arrays(nmax);
    thr[tid]->clear(nall);
  }
}

/* ---------------------------------------------------------------------- */

void FixOMP::post_force(int vflag)
{
  fprintf(stderr,"%s, vflag=%d\n", __FUNCTION__, vflag);
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(thr,force)
#endif
  {
    const int tid = get_tid();
    data_reduce_thr(&(atom->f[0][0]), nall, nthreads, 3, tid);
    if (atom->torque)
      data_reduce_thr(&(atom->torque[0][0]), nall, nthreads, 3, tid);
    //thr[tid]->reduce();
  }
}


/* ---------------------------------------------------------------------- */

double FixOMP::memory_usage()
{
  double bytes = _ndata * sizeof(ThrData *);
  bytes += _ndata * thr[0]->memory_usage();
  fprintf(stderr,"%s, bytes=%g\n", __FUNCTION__, bytes);
  
  return bytes;
}
