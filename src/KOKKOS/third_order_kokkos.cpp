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
   Contributing author: Charlie Sievers (UC Davis), charliesievers at cox.net
------------------------------------------------------------------------- */

#include "third_order_kokkos.h"

#include "angle.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "bond.h"
#include "comm.h"
#include "dihedral.h"
#include "domain.h"
#include "force.h"
#include "improper.h"
#include "kokkos.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "update.h"

using namespace LAMMPS_NS;
enum{REGULAR,ESKM};

template<class ViewA, class ViewB>
struct ForceAdder {
  ViewA a;
  ViewB b;
  ForceAdder(const ViewA& a_, const ViewB& b_):a(a_),b(b_) {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    a(i,0) += b(i,0);
    a(i,1) += b(i,1);
    a(i,2) += b(i,2);
  }
};

/* ---------------------------------------------------------------------- */

template<class View>
struct Zero {
  View v;
  Zero(const View &v_):v(v_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    v(i,0) = 0;
    v(i,1) = 0;
    v(i,2) = 0;
  }
};

/* ---------------------------------------------------------------------- */

ThirdOrderKokkos::ThirdOrderKokkos(LAMMPS *lmp) : ThirdOrder(lmp)
{
  atomKK = (AtomKokkos *) atom;
}

/* ---------------------------------------------------------------------- */

void ThirdOrderKokkos::command(int narg, char **arg)
{
  atomKK->sync(Host, X_MASK|RMASS_MASK|TYPE_MASK);
  ThirdOrder::command(narg, arg);
}

/* ----------------------------------------------------------------------
   setup without output or one-time post-init setup
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void ThirdOrderKokkos::setup()
{
  lmp->kokkos->auto_sync = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists
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
  neighbor->build(1);

  // compute all forces
  eflag=0;
  vflag=0;
  if (force->kspace) {
    force->kspace->setup();
  }
  update_force();

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);
  update->setupflag = 0;

  lmp->kokkos->auto_sync = 0;

  //if all then skip communication groupmap population
  if (gcount == atom->natoms)
    for (bigint i=0; i<atom->natoms; i++)
      groupmap[i] = i;
  else
    create_groupmap();

}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

void ThirdOrderKokkos::update_force()
{
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force_any;

  lmp->kokkos->auto_sync = 0;

  f_merge_copy = DAT::t_f_array("ThirdOrderKokkos::f_merge_copy",atomKK->k_f.extent(0));

  atomKK->modified(Host,X_MASK);
  atomKK->sync(Device,X_MASK);

  force_clear();

  neighbor->ago = 0;
  if ((modify->get_fix_by_id("package_intel")) != nullptr)
    neighbor->decide();

  if (n_pre_force) {
    modify->pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  bool execute_on_host = false;
  unsigned int datamask_read_host = 0;

  if (pair_compute_flag) {
    if (force->pair->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->pair->datamask_read;
    }
  }
  if (atomKK->molecular && force->bond)  {
    if (force->bond->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->bond->datamask_read;
    }
  }
  if (atomKK->molecular && force->angle) {
    if (force->angle->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->angle->datamask_read;
    }
  }
  if (atomKK->molecular && force->dihedral) {
    if (force->dihedral->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->dihedral->datamask_read;
    }
  }
  if (atomKK->molecular && force->improper) {
    if (force->improper->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->improper->datamask_read;
    }
  }
  if (kspace_compute_flag) {
    if (force->kspace->execution_space==Host) {
      execute_on_host  = true;
      datamask_read_host   |= force->kspace->datamask_read;
    }
  }

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    atomKK->sync(force->pair->execution_space,~(~force->pair->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    Kokkos::Timer ktimer;
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    atomKK->modified(force->pair->execution_space,~(~force->pair->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    timer->stamp(Timer::PAIR);
  }

  if (execute_on_host) {
    if (pair_compute_flag && force->pair->datamask_modify!=(F_MASK | ENERGY_MASK | VIRIAL_MASK))
      Kokkos::fence();
    atomKK->sync_overlapping_device(Host,~(~datamask_read_host|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    if (pair_compute_flag && force->pair->execution_space!=Host) {
      Kokkos::deep_copy(LMPHostType(),atomKK->k_f.h_view,0.0);
    }
  }

  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,~(~force->bond->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,~(~force->bond->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,~(~force->angle->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,~(~force->angle->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,~(~force->dihedral->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,~(~force->dihedral->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,~(~force->improper->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,~(~force->improper->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    }
    timer->stamp(Timer::BOND);
  }

  if (kspace_compute_flag) {
    atomKK->sync(force->kspace->execution_space,~(~force->kspace->datamask_read|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    force->kspace->compute(eflag,vflag);
    atomKK->modified(force->kspace->execution_space,~(~force->kspace->datamask_modify|(F_MASK | ENERGY_MASK | VIRIAL_MASK)));
    timer->stamp(Timer::KSPACE);
  }

  if (execute_on_host && !std::is_same_v<LMPHostType,LMPDeviceType>) {
    if (f_merge_copy.extent(0)<atomKK->k_f.extent(0)) {
      f_merge_copy = DAT::t_f_array("ThirdOrderKokkos::f_merge_copy",atomKK->k_f.extent(0));
    }
    f = atomKK->k_f.d_view;
    Kokkos::deep_copy(LMPHostType(),f_merge_copy,atomKK->k_f.h_view);
    Kokkos::parallel_for(atomKK->k_f.extent(0),
                         ForceAdder<DAT::t_f_array,DAT::t_f_array>(atomKK->k_f.d_view,f_merge_copy));
    atomKK->k_f.clear_sync_state(); // special case
    atomKK->k_f.modify<LMPDeviceType>();
  }
  if (n_pre_reverse) {
    modify->pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }
  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }
  // force modifications

  if (n_post_force) {
    modify->post_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  atomKK->sync(Host,F_MASK);
  lmp->kokkos->auto_sync = 1;

  ++ update->nsteps;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void ThirdOrderKokkos::force_clear()
{
  if (external_force_clear) return;

  atomKK->k_f.clear_sync_state(); // ignore host forces/torques since device views

  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  int nall = atomKK->nlocal;
  if (force->newton) nall += atomKK->nghost;

  Kokkos::parallel_for(nall, Zero<typename ArrayTypes<LMPDeviceType>::t_f_array>(atomKK->k_f.view<LMPDeviceType>()));
  atomKK->modified(Device,F_MASK);

}
