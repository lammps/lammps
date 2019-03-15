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
   OpenMP based threading support for LAMMPS
------------------------------------------------------------------------- */

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "universe.h"
#include "update.h"
#include "integrate.h"
#include "min.h"
#include "timer.h"

#include "fix_omp.h"
#include "thr_data.h"
#include "thr_omp.h"

#include "pair_hybrid.h"
#include "bond_hybrid.h"
#include "angle_hybrid.h"
#include "dihedral_hybrid.h"
#include "improper_hybrid.h"
#include "kspace.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "suffix.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static int get_tid()
{
  int tid = 0;
#if defined(_OPENMP)
  tid = omp_get_thread_num();
#endif
  return tid;
}

/* ---------------------------------------------------------------------- */

FixOMP::FixOMP(LAMMPS *lmp, int narg, char **arg)
  :  Fix(lmp, narg, arg),
     thr(NULL), last_omp_style(NULL), last_pair_hybrid(NULL),
     _nthr(-1), _neighbor(true), _mixed(false), _reduced(true)
{
  if (narg < 4) error->all(FLERR,"Illegal package omp command");

  int nthreads = 1;
  if (narg > 3) {
#if defined(_OPENMP)
    if (strcmp(arg[3],"0") == 0)
#pragma omp parallel default(none) shared(nthreads)
      nthreads = omp_get_num_threads();
    else
      nthreads = force->inumeric(FLERR,arg[3]);
#endif
  }

  if (nthreads < 1)
    error->all(FLERR,"Illegal number of OpenMP threads requested");

#if defined(_OPENMP)
  int reset_thr = 0;
#endif
  if (nthreads != comm->nthreads) {
#if defined(_OPENMP)
    reset_thr = 1;
    omp_set_num_threads(nthreads);
#endif
    comm->nthreads = nthreads;
  }

  // optional keywords

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"neigh") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package omp command");
      if (strcmp(arg[iarg+1],"yes") == 0) _neighbor = true;
      else if (strcmp(arg[iarg+1],"no") == 0) _neighbor = false;
      else error->all(FLERR,"Illegal package omp command");
      iarg += 2;
    } else error->all(FLERR,"Illegal package omp command");
  }

  // print summary of settings

  if (comm->me == 0) {
#if defined(_OPENMP)
    const char * const nmode = _neighbor ? "multi-threaded" : "serial";

    if (screen) {
      if (reset_thr)
        fprintf(screen,"set %d OpenMP thread(s) per MPI task\n", nthreads);
      fprintf(screen,"using %s neighbor list subroutines\n", nmode);
    }

    if (logfile) {
      if (reset_thr)
        fprintf(logfile,"set %d OpenMP thread(s) per MPI task\n", nthreads);
      fprintf(logfile,"using %s neighbor list subroutines\n", nmode);
    }
#else
    error->warning(FLERR,"OpenMP support not enabled during compilation; "
                         "using 1 thread only.");
#endif
  }

  // allocate list for per thread accumulator manager class instances
  // and then have each thread create an instance of this class to
  // encourage the OS to use storage that is "close" to each thread's CPU.

  thr = new ThrData *[nthreads];
  _nthr = nthreads;
#if defined(_OPENMP)
#pragma omp parallel default(none) shared(lmp)
#endif
  {
    const int tid = get_tid();
    Timer *t = new Timer(lmp);
    thr[tid] = new ThrData(tid,t);
  }
}

/* ---------------------------------------------------------------------- */

FixOMP::~FixOMP()
{
  for (int i=0; i < _nthr; ++i)
    delete thr[i];

  delete[] thr;
}

/* ---------------------------------------------------------------------- */

int FixOMP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOMP::init()
{
  // USER-OMP package cannot be used with atom_style template
  if (atom->molecular == 2)
    error->all(FLERR,"USER-OMP package does not (yet) work with "
               "atom_style template");

  // adjust number of data objects when the number of OpenMP
  // threads has been changed somehow
  const int nthreads = comm->nthreads;
  if (_nthr != nthreads) {
    if (screen) fprintf(screen,"Re-init USER-OMP for %d OpenMP thread(s)\n", nthreads);
    if (logfile) fprintf(logfile,"Re-init USER-OMP for %d OpenMP thread(s)\n", nthreads);

    for (int i=0; i < _nthr; ++i)
      delete thr[i];

    thr = new ThrData *[nthreads];
    _nthr = nthreads;
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
    {
      const int tid = get_tid();
      Timer *t = new Timer(lmp);
      thr[tid] = new ThrData(tid,t);
    }
  }

  // reset per thread timer
  for (int i=0; i < nthreads; ++i) {
    thr[i]->_timer_active=1;
    thr[i]->timer(Timer::RESET);
    thr[i]->_timer_active=-1;
  }

  if ((strstr(update->integrate_style,"respa") != NULL)
      && (strstr(update->integrate_style,"respa/omp") == NULL))
    error->all(FLERR,"Need to use respa/omp for r-RESPA with /omp styles");

  int check_hybrid, kspace_split;
  last_pair_hybrid = NULL;
  last_omp_style = NULL;
  const char *last_omp_name = NULL;
  const char *last_hybrid_name = NULL;
  const char *last_force_name = NULL;

  // support for verlet/split operation.
  // kspace_split == 0 : regular processing
  // kspace_split < 0  : master partition, does not do kspace
  // kspace_split > 0  : slave partition, only does kspace

  if (strstr(update->integrate_style,"verlet/split") != NULL) {
    if (universe->iworld == 0) kspace_split = -1;
    else kspace_split = 1;
  } else {
    kspace_split = 0;
  }

// determine which is the last force style with OpenMP
// support as this is the one that has to reduce the forces

#define CheckStyleForOMP(name)                                          \
  check_hybrid = 0;                                                     \
  if (force->name) {                                                    \
    if ( (strcmp(force->name ## _style,"hybrid") == 0) ||               \
         (strcmp(force->name ## _style,"hybrid/overlay") == 0) )        \
      check_hybrid=1;                                                   \
    if (force->name->suffix_flag & Suffix::OMP) {                       \
      last_force_name = (const char *) #name;                           \
      last_omp_name = force->name ## _style;                            \
      last_omp_style = (void *) force->name;                            \
    }                                                                   \
  }

#define CheckHybridForOMP(name,Class) \
  if (check_hybrid) {                                         \
    Class ## Hybrid *style = (Class ## Hybrid *) force->name; \
    for (int i=0; i < style->nstyles; i++) {                  \
      if (style->styles[i]->suffix_flag & Suffix::OMP) {      \
        last_force_name = (const char *) #name;               \
        last_omp_name = style->keywords[i];                   \
        last_omp_style = style->styles[i];                    \
      }                                                       \
    }                                                         \
  }

  if (kspace_split <= 0) {
    CheckStyleForOMP(pair);
    CheckHybridForOMP(pair,Pair);
    if (check_hybrid) {
      last_pair_hybrid = last_omp_style;
      last_hybrid_name = last_omp_name;
    }

    CheckStyleForOMP(bond);
    CheckHybridForOMP(bond,Bond);

    CheckStyleForOMP(angle);
    CheckHybridForOMP(angle,Angle);

    CheckStyleForOMP(dihedral);
    CheckHybridForOMP(dihedral,Dihedral);

    CheckStyleForOMP(improper);
    CheckHybridForOMP(improper,Improper);
  }

  if (kspace_split >= 0) {
    CheckStyleForOMP(kspace);
  }

#undef CheckStyleForOMP
#undef CheckHybridForOMP
  set_neighbor_omp();

  // diagnostic output
  if (comm->me == 0) {
    if (last_omp_style) {
      if (last_pair_hybrid) {
        if (screen)
          fprintf(screen,"Hybrid pair style last /omp style %s\n", last_hybrid_name);
        if (logfile)
          fprintf(logfile,"Hybrid pair style last /omp style %s\n", last_hybrid_name);
      }
      if (screen)
        fprintf(screen,"Last active /omp style is %s_style %s\n",
                last_force_name, last_omp_name);
      if (logfile)
        fprintf(logfile,"Last active /omp style is %s_style %s\n",
                last_force_name, last_omp_name);
    } else {
      if (screen)
        fprintf(screen,"No /omp style for force computation currently active\n");
      if (logfile)
        fprintf(logfile,"No /omp style for force computation currently active\n");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixOMP::set_neighbor_omp()
{
  // select or deselect multi-threaded neighbor
  // list build depending on setting in package omp.
  // NOTE: since we are at the top of the list of
  // fixes, we cannot adjust neighbor lists from
  // other fixes. those have to be re-implemented
  // as /omp fix styles. :-(

  const int neigh_omp = _neighbor ? 1 : 0;
  const int nrequest = neighbor->nrequest;

  // flag *all* neighbor list requests as USER-OMP threaded,
  // but skip lists already flagged as USER-INTEL threaded
  for (int i = 0; i < nrequest; ++i)
    if (! neighbor->requests[i]->intel)
      neighbor->requests[i]->omp = neigh_omp;
}

/* ---------------------------------------------------------------------- */

void FixOMP::setup(int)
{
  // we are post the force compute in setup. turn on timers
  for (int i=0; i < _nthr; ++i)
    thr[i]->_timer_active=0;
}

/* ---------------------------------------------------------------------- */

// adjust size and clear out per thread accumulator arrays
void FixOMP::pre_force(int)
{
  const int nall = atom->nlocal + atom->nghost;

  double **f = atom->f;
  double **torque = atom->torque;
  double *erforce = atom->erforce;
  double *de = atom->de;
  double *drho = atom->drho;

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(f,torque,erforce,de,drho)
#endif
  {
    const int tid = get_tid();
    thr[tid]->check_tid(tid);
    thr[tid]->init_force(nall,f,torque,erforce,de,drho);
  } // end of omp parallel region

  _reduced = false;
}

/* ---------------------------------------------------------------------- */

double FixOMP::memory_usage()
{
  double bytes = _nthr * (sizeof(ThrData *) + sizeof(ThrData));
  bytes += _nthr * thr[0]->memory_usage();

  return bytes;
}
