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
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "min.h"

#include "fix_omp.h"
#include "thr_data.h"
#include "thr_omp.h"

#include "pair_hybrid.h"
#include "bond_hybrid.h"
#include "angle_hybrid.h"
#include "dihedral_hybrid.h"
#include "improper_hybrid.h"
#include "kspace.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "suffix.h"

#if defined(LMP_USER_CUDA)
#include "cuda_modify_flags.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;
#if defined(LMP_USER_CUDA)
using namespace FixConstCuda;
#endif

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
              thr(NULL), last_omp_style(NULL), last_pair_hybrid(NULL),
              _nthr(-1), _neighbor(true), _newton(false)
{
  if ((narg < 4) || (narg > 6)) error->all(FLERR,"Illegal fix OMP command");
  if (strcmp(arg[1],"all") != 0) error->all(FLERR,"Illegal fix OMP command");

  int nthreads = 1;
  if (narg > 3) {
#if defined(_OPENMP)
    if (strcmp(arg[3],"*") == 0)
#pragma omp parallel default(none) shared(nthreads)
      nthreads = omp_get_num_threads();
    else
      nthreads = atoi(arg[3]);
#endif
  }

  if (nthreads < 1)
    error->all(FLERR,"Illegal number of threads requested.");

  if (nthreads != comm->nthreads) {
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
  }

  if (narg > 4) {
    if (strcmp(arg[4],"force/neigh") == 0)
      _neighbor = true;
    else if (strcmp(arg[4],"force") == 0)
      _neighbor = false;
    else
      error->all(FLERR,"Illegal fix omp mode requested.");

    if (comm->me == 0) {
      const char * const mode = _neighbor ? "OpenMP capable" : "serial";
      
      if (screen)
	fprintf(screen,"  using %s neighbor list subroutines\n", mode);
      if (logfile)
	fprintf(logfile,"  using %s neighbor list subroutines\n", mode);
    }
  }

#if 0 /* to be enabled when we can switch between half and full neighbor lists */
  if (narg > 5) {
    if (strcmp(arg[5],"neigh/half") == 0)
      _newton = true;
    else if (strcmp(arg[5],"neigh/full") == 0)
	_newton = false;
    else
      error->all(FLERR,"Illegal fix OMP command");

    if (comm->me == 0) {
      const char * const mode = _newton ? "half" : "full";
      
      if (screen)
	fprintf(screen,"  using /omp styles with %s neighbor list builds\n", mode);
      if (logfile)
	fprintf(logfile,"  using /omp styles with %s neighbor list builds\n", mode);
    }
  }
#endif

  // allocate list for per thread accumulator manager class instances
  // and then have each thread create an instance of this class to 
  // encourage the OS to use storage that is "close" to each thread's CPU.
  thr = new ThrData *[nthreads];
  _nthr = nthreads;
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
  // compatibility with USER-CUDA
  // our fix doesn't need any data transfer.
#if defined(LMP_USER_CUDA)
  if (lmp->cuda) {
    int mask = 0;
    mask |= PRE_FORCE_CUDA;
    mask |= PRE_FORCE_RESPA;
    mask |= MIN_PRE_FORCE;
    return mask;
  }
#endif

  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOMP::init()
{
  if (strstr(update->integrate_style,"respa") != NULL)
    error->all(FLERR,"Cannot use r-RESPA with /omp styles");

  int check_hybrid;
  last_pair_hybrid = NULL;
  last_omp_style = NULL;
  char *last_omp_name = NULL;
  char *last_hybrid_name = NULL;
  char *last_force_name = NULL;

// determine which is the last force style with OpenMP
// support as this is the one that has to reduce the forces

#define CheckStyleForOMP(name)						\
  check_hybrid = 0;							\
  if (force->name) {							\
    if ( (strcmp(force->name ## _style,"hybrid") == 0) ||		\
         (strcmp(force->name ## _style,"hybrid/overlay") == 0) )	\
      check_hybrid=1;							\
    if (force->name->suffix_flag & Suffix::OMP) {			\
      last_force_name = #name;						\
      last_omp_name = force->name ## _style;				\
      last_omp_style = (void *) force->name;				\
    }									\
  }

#define CheckHybridForOMP(name,Class)					\
  if (check_hybrid) {							\
    Class ## Hybrid *style = (Class ## Hybrid *) force->name;		\
    for (int i=0; i < style->nstyles; i++) {				\
      if (style->styles[i]->suffix_flag & Suffix::OMP) {		\
	last_force_name = #name;					\
	last_omp_name = style->keywords[i];				\
	last_omp_style = style->styles[i];				\
      }									\
    }									\
  }

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

  CheckStyleForOMP(kspace);

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
	fprintf(screen,"No /omp style for force computation currently active\n");
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

  for (int i = 0; i < nrequest; ++i)
    neighbor->requests[i]->omp = neigh_omp;
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
  }
}

/* ---------------------------------------------------------------------- */

double FixOMP::memory_usage()
{
  double bytes = comm->nthreads * (sizeof(ThrData *) + sizeof(ThrData));
  bytes += comm->nthreads * thr[0]->memory_usage();
  
  return bytes;
}
