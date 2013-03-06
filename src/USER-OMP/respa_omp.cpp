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
   Contributing authors: Mark Stevens (SNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "respa_omp.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
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
#include "fix_respa.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RespaOMP::RespaOMP(LAMMPS *lmp, int narg, char **arg) 
  : Respa(lmp, narg, arg),ThrOMP(lmp, THR_INTGR)
{
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void RespaOMP::init()
{
  Respa::init();

  if (atom->torque)
    error->all(FLERR,"Extended particles are not supported by respa/omp\n");
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void RespaOMP::setup()
{
  if (comm->me == 0 && screen) fprintf(screen,"Setting up run ...\n");

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
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
      {
#if defined(_OPENMP)
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	data_reduce_thr(atom->f[0], nall, nthreads, 3, tid);
      }
      fix->did_reduce();
    }
      
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  modify->setup(vflag);
  sum_flevel_f();
  output->setup();
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void RespaOMP::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear(newton[ilevel]);
    modify->setup_pre_force_respa(vflag,ilevel);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
      {
#if defined(_OPENMP)
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	data_reduce_thr(atom->f[0], nall, nthreads, 3, tid);
      }
      fix->did_reduce();
    }

    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  modify->setup(vflag);
  sum_flevel_f();
  update->setupflag = 0;
}

/* ---------------------------------------------------------------------- */

void RespaOMP::recurse(int ilevel)
{
  copy_flevel_f(ilevel);

  for (int iloop = 0; iloop < loop[ilevel]; iloop++) {

    modify->initial_integrate_respa(vflag,ilevel,iloop);
    if (modify->n_post_integrate_respa)
      modify->post_integrate_respa(ilevel,iloop);

    if (ilevel) recurse(ilevel-1);

    // at outermost level, check on rebuilding neighbor list
    // at innermost level, communicate
    // at middle levels, do nothing

    if (ilevel == nlevels-1) {
      int nflag = neighbor->decide();
      if (nflag) {
        if (modify->n_pre_exchange) modify->pre_exchange();
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
          domain->reset_box();
          comm->setup();
          if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (atom->sortfreq > 0 && 
            update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(TIME_COMM);
        if (modify->n_pre_neighbor) modify->pre_neighbor();
        neighbor->build();
        timer->stamp(TIME_NEIGHBOR);
      }

    } else if (ilevel == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(TIME_COMM);
    }

    force_clear(newton[ilevel]);
    if (modify->n_pre_force_respa)
      modify->pre_force_respa(vflag,ilevel,iloop);

    timer->stamp();
    if (level_bond == ilevel && force->bond) {
      force->bond->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_angle == ilevel && force->angle) {
      force->angle->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_dihedral == ilevel && force->dihedral) {
      force->dihedral->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_improper == ilevel && force->improper) {
      force->improper->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }
    if (level_pair == ilevel && pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }
    if (level_inner == ilevel && pair_compute_flag) {
      force->pair->compute_inner();
      timer->stamp(TIME_PAIR);
    }
    if (level_middle == ilevel && pair_compute_flag) {
      force->pair->compute_middle();
      timer->stamp(TIME_PAIR);
    }
    if (level_outer == ilevel && pair_compute_flag) {
      force->pair->compute_outer(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }
    if (level_kspace == ilevel && kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
      {
#if defined(_OPENMP)
	int tid = omp_get_thread_num();
#else
	int tid = 0;
#endif
	data_reduce_thr(atom->f[0], nall, nthreads, 3, tid);
      }
      fix->did_reduce();
    }

    if (newton[ilevel]) {
      comm->reverse_comm();
      timer->stamp(TIME_COMM);
    }

    if (modify->n_post_force_respa)
      modify->post_force_respa(vflag,ilevel,iloop);
    modify->final_integrate_respa(ilevel,iloop);
  }

  copy_f_flevel(ilevel);
}

