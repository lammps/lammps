// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mark Stevens (SNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "respa_omp.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "update.h"

#include "omp_compat.h"
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

void RespaOMP::setup(int flag)
{
  if (comm->me == 0 && screen) {
    std::string mesg = "Setting up r-RESPA/omp run ...\n";
    if (flag) {
      mesg += fmt::format("  Unit style    : {}\n", update->unit_style);
      mesg += fmt::format("  Current step  : {}\n", update->ntimestep);

      mesg += "  Time steps    :";
      for (int ilevel = 0; ilevel < nlevels; ++ilevel)
        mesg += fmt::format(" {}:{}", ilevel + 1, step[ilevel]);

      mesg += "\n  r-RESPA fixes :";
      for (int l = 0; l < modify->n_post_force_respa_any; ++l) {
        Fix *f = modify->get_fix_by_index(modify->list_post_force_respa[l]);
        if (f->respa_level >= 0)
          mesg += fmt::format(" {}:{}[{}]", MIN(f->respa_level + 1, nlevels), f->style, f->id);
      }
      mesg += "\n";
      fputs(mesg.c_str(), screen);
      timer->print_timeout(screen);
    }
  }

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
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear();
    modify->setup_pre_force_respa(vflag,ilevel);

    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
    }
    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
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

    modify->pre_reverse(eflag,vflag);
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  sum_flevel_f();
  modify->setup(vflag);
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
    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);

  for (int ilevel = 0; ilevel < nlevels; ilevel++) {
    force_clear();
    modify->setup_pre_force_respa(vflag,ilevel);

    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
    }

    if (level_pair == ilevel && pair_compute_flag)
      force->pair->compute(eflag,vflag);
    if (level_inner == ilevel && pair_compute_flag)
      force->pair->compute_inner();
    if (level_middle == ilevel && pair_compute_flag)
      force->pair->compute_middle();
    if (level_outer == ilevel && pair_compute_flag)
      force->pair->compute_outer(eflag,vflag);
    if (level_bond == ilevel && force->bond)
      force->bond->compute(eflag,vflag);
    if (level_angle == ilevel && force->angle)
      force->angle->compute(eflag,vflag);
    if (level_dihedral == ilevel && force->dihedral)
      force->dihedral->compute(eflag,vflag);
    if (level_improper == ilevel && force->improper)
      force->improper->compute(eflag,vflag);
    if (level_kspace == ilevel && force->kspace) {
      force->kspace->setup();
      if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
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

    modify->pre_reverse(eflag,vflag);
    if (newton[ilevel]) comm->reverse_comm();
    copy_f_flevel(ilevel);
  }

  sum_flevel_f();
  modify->setup(vflag);
  update->setupflag = 0;
}

/* ---------------------------------------------------------------------- */

void RespaOMP::recurse(int ilevel)
{
  copy_flevel_f(ilevel);

  for (int iloop = 0; iloop < loop[ilevel]; iloop++) {

    timer->stamp();
    modify->initial_integrate_respa(vflag,ilevel,iloop);
    if (modify->n_post_integrate_respa)
      modify->post_integrate_respa(ilevel,iloop);
    timer->stamp(Timer::MODIFY);

    // at outermost level, check on rebuilding neighbor list
    // at innermost level, communicate
    // at middle levels, do nothing

    if (ilevel == nlevels-1) {
      int nflag = neighbor->decide();
      if (nflag) {
        if (modify->n_pre_exchange) {
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
        if (atom->sortfreq > 0 &&
            update->ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (modify->n_pre_neighbor) {
          modify->pre_neighbor();
          timer->stamp(Timer::MODIFY);
        }
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
        if (modify->n_post_neighbor) {
          modify->post_neighbor();
          timer->stamp(Timer::MODIFY);
        }
      } else if (ilevel == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
      }

    } else if (ilevel == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(Timer::COMM);
    }

    // rRESPA recursion thru all levels
    // this used to be before neigh list build,
    // which prevented per-atom energy/stress being tallied correctly
    // b/c atoms migrated to new procs between short/long force calls
    // now they migrate at very start of rRESPA timestep, before all forces

    if (ilevel) recurse(ilevel-1);

    // force computations
    // important that ordering is same as Verlet
    // so that any order dependencies are the same
    // when potentials are invoked at same level

    force_clear();
    if (modify->n_pre_force_respa) {
      timer->stamp();
      modify->pre_force_respa(vflag,ilevel,iloop);
      timer->stamp(Timer::MODIFY);
    }

    timer->stamp();
    if (nhybrid_styles > 0) {
      set_compute_flags(ilevel);
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_pair == ilevel && pair_compute_flag) {
      force->pair->compute(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_inner == ilevel && pair_compute_flag) {
      force->pair->compute_inner();
      timer->stamp(Timer::PAIR);
    }
    if (level_middle == ilevel && pair_compute_flag) {
      force->pair->compute_middle();
      timer->stamp(Timer::PAIR);
    }
    if (level_outer == ilevel && pair_compute_flag) {
      force->pair->compute_outer(eflag,vflag);
      timer->stamp(Timer::PAIR);
    }
    if (level_bond == ilevel && force->bond) {
      force->bond->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_angle == ilevel && force->angle) {
      force->angle->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_dihedral == ilevel && force->dihedral) {
      force->dihedral->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_improper == ilevel && force->improper) {
      force->improper->compute(eflag,vflag);
      timer->stamp(Timer::BOND);
    }
    if (level_kspace == ilevel && kspace_compute_flag) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(Timer::KSPACE);
    }

    // reduce forces from per-thread arrays, if needed
    if (!fix->get_reduced()) {
      const int nall = atom->nlocal + atom->nghost;
      const int nthreads = comm->nthreads;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
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

    if (modify->n_pre_reverse) {
      modify->pre_reverse(eflag,vflag);
      timer->stamp(Timer::MODIFY);
    }
    if (newton[ilevel]) {
      comm->reverse_comm();
      timer->stamp(Timer::COMM);
    }
    timer->stamp();
    if (modify->n_post_force_respa_any)
      modify->post_force_respa(vflag,ilevel,iloop);
    modify->final_integrate_respa(ilevel,iloop);
    timer->stamp(Timer::MODIFY);
  }

  copy_f_flevel(ilevel);
}
