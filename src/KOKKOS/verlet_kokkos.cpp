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

#include "string.h"
#include "verlet_kokkos.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
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

#include <ctime>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletKokkos::VerletKokkos(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg) 
{
  atomKK = (AtomKokkos *) atom;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void VerletKokkos::setup()
{

  if (comm->me == 0 && screen) fprintf(screen,"Setting up run ...\n");
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atomKK->modified(Host,ALL_MASK);

  atomKK->setup();
  modify->setup_pre_exchange();
      // debug
  atomKK->sync(Host,ALL_MASK);
  atomKK->modified(Host,ALL_MASK);
  if (triclinic) domain->x2lamda(atomKK->nlocal);
  domain->pbc();

  atomKK->sync(Host,ALL_MASK);


  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();

  comm->exchange();

  if (atomKK->sortfreq > 0) atomKK->sort();

  comm->borders();

  if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);

  atomKK->sync(Host,ALL_MASK);

  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();

  atomKK->modified(Host,ALL_MASK);

  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    force->pair->compute(eflag,vflag);
    timer->stamp(TIME_PAIR);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);


  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
      force->bond->compute(eflag,vflag);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
      force->angle->compute(eflag,vflag);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
      force->dihedral->compute(eflag,vflag);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
      force->improper->compute(eflag,vflag);
    }
    timer->stamp(TIME_BOND);
  }

  if(force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    } else force->kspace->compute_dummy(eflag,vflag);
  }

  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  output->setup();
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void VerletKokkos::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    atomKK->modified(Host,ALL_MASK);

    modify->setup_pre_exchange();
      // debug
      atomKK->sync(Host,ALL_MASK);
      atomKK->modified(Host,ALL_MASK);

    if (triclinic) domain->x2lamda(atomKK->nlocal);
    domain->pbc();

    atomKK->sync(Host,ALL_MASK);

    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);

    atomKK->sync(Host,ALL_MASK);

    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();

    atomKK->modified(Host,ALL_MASK);

    neighbor->build();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    force->pair->compute(eflag,vflag);
    timer->stamp(TIME_PAIR);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);


  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
      force->bond->compute(eflag,vflag);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
      force->angle->compute(eflag,vflag);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
      force->dihedral->compute(eflag,vflag);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
      force->improper->compute(eflag,vflag);
    }
    timer->stamp(TIME_BOND);
  }

  if(force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    } else force->kspace->compute_dummy(eflag,vflag);
  }

  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void VerletKokkos::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atomKK->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  static double time = 0.0;
  static int count = 0;
  atomKK->sync(Device,ALL_MASK);
  Kokkos::Impl::Timer ktimer;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    ktimer.reset();
    modify->initial_integrate(vflag);
    time += ktimer.seconds();
    if (n_post_integrate) modify->post_integrate();

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->forward_comm();
      timer->stamp(TIME_COMM);
    } else {
      // added debug
      //atomKK->sync(Host,ALL_MASK);
      //atomKK->modified(Host,ALL_MASK);

      if (n_pre_exchange) modify->pre_exchange();
      // debug
      //atomKK->sync(Host,ALL_MASK);
      //atomKK->modified(Host,ALL_MASK);
      if (triclinic) domain->x2lamda(atomKK->nlocal);
      domain->pbc();
      if (domain->box_change) {
        domain->reset_box();
        comm->setup();
        if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();

      // added debug
      //atomKK->sync(Device,ALL_MASK);
      //atomKK->modified(Device,ALL_MASK);

      comm->exchange();
      if (sortflag && ntimestep >= atomKK->nextsort) atomKK->sort();
      comm->borders();

      // added debug
      //atomKK->sync(Host,ALL_MASK);
      //atomKK->modified(Host,ALL_MASK);

      if (triclinic) domain->lamda2x(atomKK->nlocal+atomKK->nghost);

      timer->stamp(TIME_COMM);
      if (n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
      timer->stamp(TIME_NEIGHBOR);
    }

    // force computations
    // important for pair to come before bonded contributions
    // since some bonded potentials tally pairwise energy/virial
    // and Pair:ev_tally() needs to be called before any tallying

    force_clear();
    // added for debug
    //atomKK->k_x.sync<LMPHostType>();
    //atomKK->k_f.sync<LMPHostType>();
    //atomKK->k_f.modify<LMPHostType>();
    if (n_pre_force) modify->pre_force(vflag);

    timer->stamp();

    if (pair_compute_flag) {
      atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
      atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }

    if (atomKK->molecular) {
      if (force->bond) {
        atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
        atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
        force->bond->compute(eflag,vflag);
      }
      if (force->angle) {
        atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
        atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
        force->angle->compute(eflag,vflag);
      }
      if (force->dihedral) {
        atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
        atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
        force->dihedral->compute(eflag,vflag);
      }
      if (force->improper) {
        atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
        atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
        force->improper->compute(eflag,vflag);
      }
      timer->stamp(TIME_BOND);
    }

    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    // reverse communication of forces

    if (force->newton) comm->reverse_comm();
    timer->stamp(TIME_COMM);

    // force modifications, final time integration, diagnostics

    ktimer.reset();

    if (n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (n_end_of_step) modify->end_of_step();

    time += ktimer.seconds();

    // all output

    if (ntimestep == output->next) {
       atomKK->sync(Host,ALL_MASK);

      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void VerletKokkos::force_clear()
{
  int i;

  if (external_force_clear) return;

  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  if (neighbor->includegroup == 0) {
    int nall;
    if (force->newton) nall = atomKK->nlocal + atomKK->nghost;
    else nall = atomKK->nlocal;

    size_t nbytes = sizeof(double) * nall;

    if (nbytes) {
      if (atomKK->k_f.modified_host > atomKK->k_f.modified_device) {
    	memset_kokkos(atomKK->k_f.view<LMPHostType>());
    	atomKK->modified(Host,F_MASK);
      } else {
        memset_kokkos(atomKK->k_f.view<LMPDeviceType>());
        atomKK->modified(Device,F_MASK);
      }
      if (torqueflag)  memset(&(atomKK->torque[0][0]),0,3*nbytes);

    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    int nall = atomKK->nfirst;
    if (atomKK->k_f.modified_host > atomKK->k_f.modified_device) {
      memset_kokkos(atomKK->k_f.view<LMPHostType>());
      atomKK->modified(Host,F_MASK);
    } else {
      memset_kokkos(atomKK->k_f.view<LMPDeviceType>());
      atomKK->modified(Device,F_MASK);
    }
    if (torqueflag) {
      double **torque = atomKK->torque;
      for (i = 0; i < nall; i++) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
        torque[i][2] = 0.0;
      }
    }

    if (force->newton) {
      nall = atomKK->nlocal + atomKK->nghost;

      if (torqueflag) {
        double **torque = atomKK->torque;
        for (i = atomKK->nlocal; i < nall; i++) {
          torque[i][0] = 0.0;
          torque[i][1] = 0.0;
          torque[i][2] = 0.0;
        }
      }

    }
  }
}
