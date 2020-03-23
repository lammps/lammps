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

#include <cstring>
#include "verlet_lrt_intel.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletLRTIntel::VerletLRTIntel(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg) {
  #if defined(_LMP_INTEL_LRT_PTHREAD)
  pthread_mutex_init(&_kmutex,NULL);
  #endif
}

/* ---------------------------------------------------------------------- */

VerletLRTIntel::~VerletLRTIntel()
{
  #if defined(_LMP_INTEL_LRT_PTHREAD)
  pthread_mutex_destroy(&_kmutex);
  #endif
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletLRTIntel::init()
{
  Verlet::init();

  _intel_kspace = (PPPMIntel*)(force->kspace_match("^pppm/intel", 0));

  #ifndef LMP_INTEL_USELRT
  error->all(FLERR,
             "LRT otion for Intel package disabled at compile time");
  #endif
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void VerletLRTIntel::setup(int flag)
{
  if (_intel_kspace == 0) {
    Verlet::setup(flag);
    return;
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_intel_kspace->use_base()) {
    _intel_kspace = 0;
    Verlet::setup(flag);
    return;
  }
  #endif

  if (comm->me == 0 && screen) {
    fprintf(screen,"Setting up VerletLRTIntel run ...\n");
    fprintf(screen,"  Unit style    : %s\n", update->unit_style);
    fprintf(screen,"  Current step  : " BIGINT_FORMAT "\n", update->ntimestep);
    fprintf(screen,"  Time step     : %g\n", update->dt);
    timer->print_timeout(screen);
  }

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  #if defined(__linux)
  if (comm->me == 0) {
    cpu_set_t cpuset;
    sched_getaffinity(0, sizeof(cpuset), &cpuset);
    int my_cpu_count = CPU_COUNT(&cpuset);
    if (my_cpu_count < comm->nthreads + 1) {
      char str[128];
      sprintf(str,"Using %d threads per MPI rank, but only %d core(s)"
                  " allocated for each MPI rank",
              comm->nthreads + 1, my_cpu_count);
      error->warning(FLERR, str);
    }
  }
  #endif

  _kspace_ready = 0;
  _kspace_done = 0;
  pthread_cond_init(&_kcond, NULL);
  pthread_attr_init(&_kspace_attr);
  pthread_attr_setdetachstate(&_kspace_attr, PTHREAD_CREATE_JOINABLE);
  #endif

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  modify->setup_pre_exchange();
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  neighbor->build(1);
  neighbor->ncalls = 0;

  // compute all forces

  force->setup();
  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);
  _intel_kspace->setup();

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  pthread_create(&_kspace_thread, &_kspace_attr,
                 &VerletLRTIntel::k_launch_loop, this);
  #elif defined(_LMP_INTEL_LRT_11)
  std::thread _kspace_thread;
  if (kspace_compute_flag)
    _kspace_thread=std::thread([=]{ _intel_kspace->compute_first(eflag,
                                                                 vflag); });
  else
    _kspace_thread=std::thread([=]{ _intel_kspace->compute_dummy(eflag,
                                                                 vflag); });
  #endif

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  pthread_mutex_lock(&_kmutex);
  while (!_kspace_done)
    pthread_cond_wait(&_kcond, &_kmutex);
  _kspace_done = 0;
  pthread_mutex_unlock(&_kmutex);
  #elif defined(_LMP_INTEL_LRT_11)
  _kspace_thread.join();
  #endif

  if (kspace_compute_flag) _intel_kspace->compute_second(eflag,vflag);

  modify->pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  output->setup();
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void VerletLRTIntel::run(int n)
{
  if (_intel_kspace == 0) {
    Verlet::run(n);
    return;
  }

  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  _krun_n = n;
  #endif

  int run_cancelled = 0;

  for (int i = 0; i < n; i++) {
    if (timer->check_timeout(i)) {
      update->nsteps = i;
      run_cancelled = 1;
      break;
    }

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    timer->stamp();
    modify->initial_integrate(vflag);
    if (n_post_integrate) modify->post_integrate();
    timer->stamp(Timer::MODIFY);

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(Timer::COMM);
      _intel_kspace->pack_buffers();
    } else {
      if (n_pre_exchange) {
        timer->stamp();
        modify->pre_exchange();
        timer->stamp(Timer::MODIFY);
      }
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      if (domain->box_change) {
        domain->reset_box();
        comm->setup();
        if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();
      comm->exchange();
      if (sortflag && ntimestep >= atom->nextsort) atom->sort();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      timer->stamp(Timer::COMM);
      if (n_pre_neighbor) {
        modify->pre_neighbor();
        timer->stamp(Timer::MODIFY);
      }
      neighbor->build(1);
      timer->stamp(Timer::NEIGH);
    }

    // force computations
    // important for pair to come before bonded contributions
    // since some bonded potentials tally pairwise energy/virial
    // and Pair:ev_tally() needs to be called before any tallying

    force_clear();

    timer->stamp();

    #if defined(_LMP_INTEL_LRT_PTHREAD)
    pthread_mutex_lock(&_kmutex);
    _kspace_ready = 1;
    _kspace_done = 0;
    pthread_cond_signal(&_kcond);
    pthread_mutex_unlock(&_kmutex);
    #elif defined(_LMP_INTEL_LRT_11)
    std::thread _kspace_thread;
    if (kspace_compute_flag)
      _kspace_thread=std::thread([=] {
        _intel_kspace->compute_first(eflag, vflag);
        timer->stamp(Timer::KSPACE);
      } );
    #endif

    if (n_pre_force) {
      modify->pre_force(vflag);
      timer->stamp(Timer::MODIFY);
    }

    if (pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }

    if (atom->molecular) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }

    #if defined(_LMP_INTEL_LRT_PTHREAD)
    pthread_mutex_lock(&_kmutex);
    while (!_kspace_done)
      pthread_cond_wait(&_kcond, &_kmutex);
    _kspace_done = 0;
    pthread_mutex_unlock(&_kmutex);
    #elif defined(_LMP_INTEL_LRT_11)
    if (kspace_compute_flag)
      _kspace_thread.join();
    #endif

    if (kspace_compute_flag) {
      _intel_kspace->compute_second(eflag,vflag);
      timer->stamp(Timer::KSPACE);
    }

    if (n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);
    }

    // reverse communication of forces

    if (force->newton) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
    }

    // force modifications, final time integration, diagnostics

    if (n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (n_end_of_step) modify->end_of_step();
    timer->stamp(Timer::MODIFY);

    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  #if defined(_LMP_INTEL_LRT_PTHREAD)
  if (run_cancelled)
    pthread_cancel(_kspace_thread);
  else {
    pthread_mutex_lock(&_kmutex);
    _kspace_ready = 1;
    _kspace_done = 0;
    pthread_cond_signal(&_kcond);
    pthread_mutex_unlock(&_kmutex);
    pthread_join(_kspace_thread, NULL);
    pthread_attr_destroy(&_kspace_attr);
  }
  #endif
}

#if defined(_LMP_INTEL_LRT_PTHREAD)

/* ----------------------------------------------------------------------
   PTHREAD Loop for long-range
------------------------------------------------------------------------- */
void * VerletLRTIntel::k_launch_loop(void *context)
{
  VerletLRTIntel * const c = (VerletLRTIntel *)context;

  if (c->kspace_compute_flag)
    c->_intel_kspace->compute_first(c->eflag, c->vflag);
  else
    c->_intel_kspace->compute_dummy(c->eflag, c->vflag);

  pthread_mutex_lock(&(c->_kmutex));
  c->_kspace_done = 1;
  pthread_cond_signal(&(c->_kcond));
  pthread_mutex_unlock(&(c->_kmutex));

  pthread_mutex_lock(&(c->_kmutex));
  while (!(c->_kspace_ready))
    pthread_cond_wait(&(c->_kcond), &(c->_kmutex));
  c->_kspace_ready = 0;
  const int n = c->_krun_n;
  pthread_mutex_unlock(&(c->_kmutex));

  for (int i = 0; i < n; i++) {

    if (c->kspace_compute_flag) {
      c->_intel_kspace->compute_first(c->eflag, c->vflag);
      c->timer->stamp(Timer::KSPACE);
    }

    pthread_mutex_lock(&(c->_kmutex));
    c->_kspace_done = 1;
    pthread_cond_signal(&(c->_kcond));
    pthread_mutex_unlock(&(c->_kmutex));

    pthread_mutex_lock(&(c->_kmutex));
    while (!(c->_kspace_ready))
      pthread_cond_wait(&(c->_kcond), &(c->_kmutex));
    c->_kspace_ready = 0;
    pthread_mutex_unlock(&(c->_kmutex));
  }

  pthread_exit(NULL);
  return NULL;
}

#endif
