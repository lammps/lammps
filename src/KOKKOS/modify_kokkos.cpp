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

#include "modify_kokkos.h"
#include "atom_kokkos.h"
#include "update.h"
#include "fix.h"
#include "compute.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ModifyKokkos::ModifyKokkos(LAMMPS *lmp) : Modify(lmp)
{
  atomKK = (AtomKokkos *) atom;
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup(int vflag)
{
  // compute setup needs to come before fix setup
  //   b/c NH fixes need DOF of temperature computes
  // fix group setup() is special case since populates a dynamic group
  //   needs to be done before temperature compute setup

  for (int i = 0; i < nfix; i++) {
    if (strcmp(fix[i]->style,"GROUP") == 0) {
      atomKK->sync(fix[i]->execution_space,fix[i]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[i]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[i]->setup(vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[i]->execution_space,fix[i]->datamask_modify);
    }
  }

  for (int i = 0; i < ncompute; i++) compute[i]->setup();

  if (update->whichflag == 1)
    for (int i = 0; i < nfix; i++) {
      atomKK->sync(fix[i]->execution_space,fix[i]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[i]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[i]->setup(vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[i]->execution_space,fix[i]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < nfix; i++) {
      atomKK->sync(fix[i]->execution_space,fix[i]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[i]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[i]->min_setup(vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[i]->execution_space,fix[i]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    for (int i = 0; i < n_pre_exchange; i++) {
      atomKK->sync(fix[list_pre_exchange[i]]->execution_space,
                   fix[list_pre_exchange[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_pre_exchange[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_pre_exchange[i]]->setup_pre_exchange();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_pre_exchange[i]]->execution_space,
                       fix[list_pre_exchange[i]]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_exchange; i++) {
      atomKK->sync(fix[list_min_pre_exchange[i]]->execution_space,
                   fix[list_min_pre_exchange[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_min_pre_exchange[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_min_pre_exchange[i]]->setup_pre_exchange();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_min_pre_exchange[i]]->execution_space,
                       fix[list_min_pre_exchange[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_neighbor; i++) {
      atomKK->sync(fix[list_pre_neighbor[i]]->execution_space,
                   fix[list_pre_neighbor[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_pre_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_pre_neighbor[i]]->setup_pre_neighbor();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_pre_neighbor[i]]->execution_space,
                       fix[list_pre_neighbor[i]]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_neighbor; i++) {
      atomKK->sync(fix[list_min_pre_neighbor[i]]->execution_space,
                   fix[list_min_pre_neighbor[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_min_pre_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_min_pre_neighbor[i]]->setup_pre_neighbor();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_min_pre_neighbor[i]]->execution_space,
                       fix[list_min_pre_neighbor[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   setup post_neighbor call, only for fixes that define post_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyKokkos::setup_post_neighbor()
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_post_neighbor; i++) {
      atomKK->sync(fix[list_post_neighbor[i]]->execution_space,
                   fix[list_post_neighbor[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_post_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_post_neighbor[i]]->setup_post_neighbor();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_post_neighbor[i]]->execution_space,
                       fix[list_post_neighbor[i]]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_post_neighbor; i++) {
      atomKK->sync(fix[list_min_post_neighbor[i]]->execution_space,
                   fix[list_min_post_neighbor[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_min_post_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_min_post_neighbor[i]]->setup_post_neighbor();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_min_post_neighbor[i]]->execution_space,
                       fix[list_min_post_neighbor[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force(int vflag)
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_force; i++) {
      atomKK->sync(fix[list_pre_force[i]]->execution_space,
                   fix[list_pre_force[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_pre_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_pre_force[i]]->setup_pre_force(vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_pre_force[i]]->execution_space,
                       fix[list_pre_force[i]]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_force; i++) {
      atomKK->sync(fix[list_min_pre_force[i]]->execution_space,
                   fix[list_min_pre_force[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_min_pre_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_min_pre_force[i]]->setup_pre_force(vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_min_pre_force[i]]->execution_space,
                       fix[list_min_pre_force[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_reverse(int eflag, int vflag)
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_reverse; i++) {
      atomKK->sync(fix[list_pre_reverse[i]]->execution_space,
                   fix[list_pre_reverse[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_pre_reverse[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_pre_reverse[i]]->setup_pre_reverse(eflag,vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_pre_reverse[i]]->execution_space,
                       fix[list_pre_reverse[i]]->datamask_modify);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_reverse; i++) {
      atomKK->sync(fix[list_min_pre_reverse[i]]->execution_space,
                   fix[list_min_pre_reverse[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_min_pre_reverse[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_min_pre_reverse[i]]->setup_pre_reverse(eflag,vflag);
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_min_pre_reverse[i]]->execution_space,
                       fix[list_min_pre_reverse[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate(int vflag)
{
  for (int i = 0; i < n_initial_integrate; i++) {
    atomKK->sync(fix[list_initial_integrate[i]]->execution_space,
                 fix[list_initial_integrate[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_initial_integrate[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_initial_integrate[i]]->initial_integrate(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_initial_integrate[i]]->execution_space,
                     fix[list_initial_integrate[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate()
{
  for (int i = 0; i < n_post_integrate; i++) {
    atomKK->sync(fix[list_post_integrate[i]]->execution_space,
                 fix[list_post_integrate[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_post_integrate[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_post_integrate[i]]->post_integrate();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_post_integrate[i]]->execution_space,
                     fix[list_post_integrate[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_exchange()
{
  for (int i = 0; i < n_pre_exchange; i++) {
    atomKK->sync(fix[list_pre_exchange[i]]->execution_space,
                 fix[list_pre_exchange[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_exchange[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_exchange[i]]->pre_exchange();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_exchange[i]]->execution_space,
                     fix[list_pre_exchange[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_neighbor()
{
  for (int i = 0; i < n_pre_neighbor; i++) {
    atomKK->sync(fix[list_pre_neighbor[i]]->execution_space,
                 fix[list_pre_neighbor[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_neighbor[i]]->pre_neighbor();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_neighbor[i]]->execution_space,
                     fix[list_pre_neighbor[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   post_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_neighbor()
{
  for (int i = 0; i < n_post_neighbor; i++) {
    atomKK->sync(fix[list_post_neighbor[i]]->execution_space,
                 fix[list_post_neighbor[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_post_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_post_neighbor[i]]->post_neighbor();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_post_neighbor[i]]->execution_space,
                     fix[list_post_neighbor[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force(int vflag)
{
  for (int i = 0; i < n_pre_force; i++) {
    atomKK->sync(fix[list_pre_force[i]]->execution_space,
                 fix[list_pre_force[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_force[i]]->pre_force(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_force[i]]->execution_space,
                     fix[list_pre_force[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_reverse(int eflag, int vflag)
{
  for (int i = 0; i < n_pre_reverse; i++) {
    atomKK->sync(fix[list_pre_reverse[i]]->execution_space,
                 fix[list_pre_reverse[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_reverse[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_reverse[i]]->pre_reverse(eflag,vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_reverse[i]]->execution_space,
                     fix[list_pre_reverse[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_force(int vflag)
{
  for (int i = 0; i < n_post_force; i++) {
    atomKK->sync(fix[list_post_force[i]]->execution_space,
                 fix[list_post_force[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_post_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_post_force[i]]->post_force(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_post_force[i]]->execution_space,
                     fix[list_post_force[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++) {
    atomKK->sync(fix[list_final_integrate[i]]->execution_space,
                 fix[list_final_integrate[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_final_integrate[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_final_integrate[i]]->final_integrate();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_final_integrate[i]]->execution_space,
                     fix[list_final_integrate[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   fused initial and final integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::fused_integrate(int vflag)
{
  for (int i = 0; i < n_final_integrate; i++) {
    atomKK->sync(fix[list_final_integrate[i]]->execution_space,
                 fix[list_final_integrate[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_final_integrate[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_final_integrate[i]]->fused_integrate(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_final_integrate[i]]->execution_space,
                     fix[list_final_integrate[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyKokkos::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0) {
      atomKK->sync(fix[list_end_of_step[i]]->execution_space,
                   fix[list_end_of_step[i]]->datamask_read);
      int prev_auto_sync = lmp->kokkos->auto_sync;
      if (!fix[list_end_of_step[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
      fix[list_end_of_step[i]]->end_of_step();
      lmp->kokkos->auto_sync = prev_auto_sync;
      atomKK->modified(fix[list_end_of_step[i]]->execution_space,
                       fix[list_end_of_step[i]]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   coupling energy call, only for relevant fixes
   each thermostsat fix returns this via compute_scalar()
   ecouple = cumulative energy added to reservoir by thermostatting
------------------------------------------------------------------------- */

double ModifyKokkos::energy_couple()
{
  double energy = 0.0;
  for (int i = 0; i < n_energy_couple; i++) {
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_energy_couple[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    energy += fix[list_energy_couple[i]]->compute_scalar();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_energy_couple[i]]->execution_space,
                     fix[list_energy_couple[i]]->datamask_modify);
  }
  return energy;
}

/* ----------------------------------------------------------------------
   global energy call, only for relevant fixes
   they return energy via compute_scalar()
   called by compute pe
------------------------------------------------------------------------- */

double ModifyKokkos::energy_global()
{
  double energy = 0.0;
  for (int i = 0; i < n_energy_global; i++) {
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_energy_global[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    energy += fix[list_energy_global[i]]->compute_scalar();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_energy_global[i]]->execution_space,
                     fix[list_energy_global[i]]->datamask_modify);
  }
  return energy;
}

/* ----------------------------------------------------------------------
   peratom energy call, only for relevant fixes
   called by compute pe/atom
------------------------------------------------------------------------- */

void ModifyKokkos::energy_atom(int nlocal, double *energy)
{
  int i,j;
  double *eatom;

  for (i = 0; i < n_energy_atom; i++) {
    eatom = fix[list_energy_atom[i]]->eatom;
    if (!eatom) continue;
    for (j = 0; j < nlocal; j++) energy[j] += eatom[j];
  }
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void ModifyKokkos::post_run()
{
  for (int i = 0; i < nfix; i++) {
    atomKK->sync(fix[i]->execution_space,
                 fix[i]->datamask_read);
    fix[i]->post_run();
    atomKK->modified(fix[i]->execution_space,
                     fix[i]->datamask_modify);
  }

  // must reset this to its default value, since computes may be added
  // or removed between runs and with this change we will redirect any
  // calls to addstep_compute() to addstep_compute_all() instead.
  n_timeflag = -1;

}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force_respa(int vflag, int ilevel)
{
  for (int i = 0; i < n_pre_force; i++) {
    atomKK->sync(fix[list_pre_force[i]]->execution_space,
                 fix[list_pre_force[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_force[i]]->setup_pre_force_respa(vflag,ilevel);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_force[i]]->execution_space,
                     fix[list_pre_force[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_initial_integrate_respa; i++) {
    atomKK->sync(fix[list_initial_integrate_respa[i]]->execution_space,
                 fix[list_initial_integrate_respa[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_initial_integrate_respa[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_initial_integrate_respa[i]]->
      initial_integrate_respa(vflag,ilevel,iloop);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_initial_integrate_respa[i]]->execution_space,
                     fix[list_initial_integrate_respa[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_post_integrate_respa; i++) {
    atomKK->sync(fix[list_post_integrate_respa[i]]->execution_space,
                 fix[list_post_integrate_respa[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_post_integrate_respa[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_post_integrate_respa[i]]->post_integrate_respa(ilevel,iloop);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_post_integrate_respa[i]]->execution_space,
                     fix[list_post_integrate_respa[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_pre_force_respa; i++) {
    atomKK->sync(fix[list_pre_force_respa[i]]->execution_space,
                 fix[list_pre_force_respa[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_pre_force_respa[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_pre_force_respa[i]]->pre_force_respa(vflag,ilevel,iloop);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_pre_force_respa[i]]->execution_space,
                     fix[list_pre_force_respa[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_post_force_respa; i++) {
    atomKK->sync(fix[list_post_force_respa[i]]->execution_space,
                 fix[list_post_force_respa[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_post_force_respa[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_post_force_respa[i]]->post_force_respa(vflag,ilevel,iloop);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_post_force_respa[i]]->execution_space,
                     fix[list_post_force_respa[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_final_integrate_respa; i++) {
    atomKK->sync(fix[list_final_integrate_respa[i]]->execution_space,
                 fix[list_final_integrate_respa[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_final_integrate_respa[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_final_integrate_respa[i]]->final_integrate_respa(ilevel,iloop);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_final_integrate_respa[i]]->execution_space,
                     fix[list_final_integrate_respa[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_exchange()
{
  for (int i = 0; i < n_min_pre_exchange; i++) {
    atomKK->sync(fix[list_min_pre_exchange[i]]->execution_space,
                 fix[list_min_pre_exchange[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_pre_exchange[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_pre_exchange[i]]->min_pre_exchange();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_pre_exchange[i]]->execution_space,
                     fix[list_min_pre_exchange[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_neighbor()
{
  for (int i = 0; i < n_min_pre_neighbor; i++) {
    atomKK->sync(fix[list_min_pre_neighbor[i]]->execution_space,
                 fix[list_min_pre_neighbor[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_pre_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_pre_neighbor[i]]->min_pre_neighbor();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_pre_neighbor[i]]->execution_space,
                     fix[list_min_pre_neighbor[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer post-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_post_neighbor()
{
  for (int i = 0; i < n_min_post_neighbor; i++) {
    atomKK->sync(fix[list_min_post_neighbor[i]]->execution_space,
                 fix[list_min_post_neighbor[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_post_neighbor[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_post_neighbor[i]]->min_post_neighbor();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_post_neighbor[i]]->execution_space,
                     fix[list_min_post_neighbor[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_force(int vflag)
{
  for (int i = 0; i < n_min_pre_force; i++) {
    atomKK->sync(fix[list_min_pre_force[i]]->execution_space,
                 fix[list_min_pre_force[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_pre_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_pre_force[i]]->min_pre_force(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_pre_force[i]]->execution_space,
                     fix[list_min_pre_force[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_reverse(int eflag, int vflag)
{
  for (int i = 0; i < n_min_pre_reverse; i++) {
    atomKK->sync(fix[list_min_pre_reverse[i]]->execution_space,
                 fix[list_min_pre_reverse[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_pre_reverse[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_pre_reverse[i]]->min_pre_reverse(eflag,vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_pre_reverse[i]]->execution_space,
                     fix[list_min_pre_reverse[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_post_force(int vflag)
{
  for (int i = 0; i < n_min_post_force; i++) {
    atomKK->sync(fix[list_min_post_force[i]]->execution_space,
                 fix[list_min_post_force[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_post_force[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_post_force[i]]->min_post_force(vflag);
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_post_force[i]]->execution_space,
                     fix[list_min_post_force[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double ModifyKokkos::min_energy(double *fextra)
{
  int ifix,index;

  index = 0;
  double eng = 0.0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[ifix]->kokkosable) lmp->kokkos->auto_sync = 1;
    eng += fix[ifix]->min_energy(&fextra[index]);
    index += fix[ifix]->min_dof();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
  }
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_store()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_energy[i]]->min_store();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   manage state of extra minimizer dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_clearstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_energy[i]]->min_clearstore();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
}

void ModifyKokkos::min_pushstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_energy[i]]->min_pushstore();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
}

void ModifyKokkos::min_popstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[list_min_energy[i]]->min_popstore();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   displace extra minimizer dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_step(double alpha, double *hextra)
{
  int ifix,index;

  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[ifix]->kokkosable) lmp->kokkos->auto_sync = 1;
    fix[ifix]->min_step(alpha,&hextra[index]);
    index += fix[ifix]->min_dof();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double ModifyKokkos::max_alpha(double *hextra)
{
  int ifix,index;

  double alpha = BIG;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[ifix]->kokkosable) lmp->kokkos->auto_sync = 1;
    double alpha_one = fix[ifix]->max_alpha(&hextra[index]);
    alpha = MIN(alpha,alpha_one);
    index += fix[ifix]->min_dof();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_dof()
{
  int ndof = 0;
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    ndof += fix[list_min_energy[i]]->min_dof();
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
  return ndof;
}

/* ----------------------------------------------------------------------
   reset minimizer reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_reset_ref()
{
  int itmp,itmpall;
  itmpall = 0;
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    if (!fix[list_min_energy[i]]->kokkosable) lmp->kokkos->auto_sync = 1;
    itmp = fix[list_min_energy[i]]->min_reset_ref();
    if (itmp) itmpall = 1;
    lmp->kokkos->auto_sync = prev_auto_sync;
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
  }
  return itmpall;
}

/* ----------------------------------------------------------------------
   check if initial and final integrate can be fused
------------------------------------------------------------------------- */

int ModifyKokkos::check_fuse_integrate()
{
  int fuse_integrate_flag = 1;

  for (int i = 0; i < n_initial_integrate; i++)
    if (!fix[list_initial_integrate[i]]->fuse_integrate_flag)
      fuse_integrate_flag = 0;

  for (int i = 0; i < n_final_integrate; i++)
    if (!fix[list_final_integrate[i]]->fuse_integrate_flag)
      fuse_integrate_flag = 0;

  return fuse_integrate_flag;
}
