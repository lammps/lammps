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
#include "compute.h"
#include "fix.h"
#include "kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

ModifyKokkos::ModifyKokkos(LAMMPS *lmp) : Modify(lmp)
{
  atomKK = (AtomKokkos *) atom;
  fuse_integrate_ok = 0;
}

/* ----------------------------------------------------------------------
   invoke func(obj) with per-atom data synchronized to the execution space
   of obj, and with automatic synchronization on for non-KOKKOS styles
------------------------------------------------------------------------- */

template <typename T, typename Func> void ModifyKokkos::sync_call(T *obj, Func &&func)
{
  atomKK->sync(obj->execution_space, obj->datamask_read);
  int prev_auto_sync = lmp->kokkos->auto_sync;
  if (!obj->kokkosable) lmp->kokkos->auto_sync = 1;
  func(obj);
  lmp->kokkos->auto_sync = prev_auto_sync;
  atomKK->modified(obj->execution_space, obj->datamask_modify);
}

template <typename Func> void ModifyKokkos::for_each_fix_kk(int bit, Func &&func)
{
  for (std::size_t i = 0; i < fixes.size(); ++i)
    if (callmask[i] & bit) sync_call(fixes[i], func);
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes, then cache whether the initial and
   final integrate calls can be fused: this is a property of the fixes which
   does not change between inits, so it is not re-checked on every step
------------------------------------------------------------------------- */

void ModifyKokkos::init()
{
  Modify::init();

  fuse_integrate_ok = 1;
  for (std::size_t i = 0; i < fixes.size(); ++i)
    if ((callmask[i] & (INITIAL_INTEGRATE | FINAL_INTEGRATE)) && !fixes[i]->fuse_integrate_flag)
      fuse_integrate_ok = 0;
}

/* ----------------------------------------------------------------------
   check if initial and final integrate can be fused
------------------------------------------------------------------------- */

int ModifyKokkos::check_fuse_integrate()
{
  return fuse_integrate_ok;
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

  for_each_fix_kk(GROUP_FIX, [vflag](Fix *ifix) {
    ifix->setup(vflag);
  });

  // the computes get the same treatment as the fixes: several of them
  // read the per-atom arrays in setup(), through the plain pointers when they
  // have no KOKKOS version, and were reaching them without a transfer

  for (std::size_t i = 0; i < computes.size(); i++)
    sync_call(computes[i], [](Compute *icompute) {
      icompute->setup();
    });

  if (update->whichflag == 1)
    for (std::size_t i = 0; i < fixes.size(); i++)
      sync_call(fixes[i], [vflag](Fix *ifix) {
        ifix->setup(vflag);
      });
  else if (update->whichflag == 2)
    for (std::size_t i = 0; i < fixes.size(); i++)
      sync_call(fixes[i], [vflag](Fix *ifix) {
        ifix->min_setup(vflag);
      });
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    for_each_fix_kk(PRE_EXCHANGE, [](Fix *ifix) {
      ifix->setup_pre_exchange();
    });
  else if (update->whichflag == 2)
    for_each_fix_kk(MIN_PRE_EXCHANGE, [](Fix *ifix) {
      ifix->setup_pre_exchange();
    });
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    for_each_fix_kk(PRE_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_pre_neighbor();
    });
  else if (update->whichflag == 2)
    for_each_fix_kk(MIN_PRE_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_pre_neighbor();
    });
}

/* ----------------------------------------------------------------------
   setup post_neighbor call, only for fixes that define post_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyKokkos::setup_post_neighbor()
{
  if (update->whichflag == 1)
    for_each_fix_kk(POST_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_post_neighbor();
    });
  else if (update->whichflag == 2)
    for_each_fix_kk(MIN_POST_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_post_neighbor();
    });
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force(int vflag)
{
  if (update->whichflag == 1)
    for_each_fix_kk(PRE_FORCE, [vflag](Fix *ifix) {
      ifix->setup_pre_force(vflag);
    });
  else if (update->whichflag == 2)
    for_each_fix_kk(MIN_PRE_FORCE, [vflag](Fix *ifix) {
      ifix->setup_pre_force(vflag);
    });
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_reverse(int eflag, int vflag)
{
  if (update->whichflag == 1)
    for_each_fix_kk(PRE_REVERSE, [=](Fix *ifix) {
      ifix->setup_pre_reverse(eflag, vflag);
    });
  else if (update->whichflag == 2)
    for_each_fix_kk(MIN_PRE_REVERSE, [=](Fix *ifix) {
      ifix->setup_pre_reverse(eflag, vflag);
    });
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate(int vflag)
{
  for_each_fix_kk(INITIAL_INTEGRATE, [vflag](Fix *ifix) {
    ifix->initial_integrate(vflag);
  });
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate()
{
  for_each_fix_kk(POST_INTEGRATE, [](Fix *ifix) {
    ifix->post_integrate();
  });
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_exchange()
{
  for_each_fix_kk(PRE_EXCHANGE, [](Fix *ifix) {
    ifix->pre_exchange();
  });
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_neighbor()
{
  for_each_fix_kk(PRE_NEIGHBOR, [](Fix *ifix) {
    ifix->pre_neighbor();
  });
}

/* ----------------------------------------------------------------------
   post_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_neighbor()
{
  for_each_fix_kk(POST_NEIGHBOR, [](Fix *ifix) {
    ifix->post_neighbor();
  });
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force(int vflag)
{
  for_each_fix_kk(PRE_FORCE, [vflag](Fix *ifix) {
    ifix->pre_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_reverse(int eflag, int vflag)
{
  for_each_fix_kk(PRE_REVERSE, [=](Fix *ifix) {
    ifix->pre_reverse(eflag, vflag);
  });
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
------------------------------------------------------------------------- */

void ModifyKokkos::post_force(int vflag)
{
  if (n_post_force_group)
    for_each_fix_kk(GROUP_FIX, [vflag](Fix *ifix) {
      ifix->post_force(vflag);
    });
  if (n_post_force)
    for_each_fix_kk(POST_FORCE, [vflag](Fix *ifix) {
      ifix->post_force(vflag);
    });
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate()
{
  for_each_fix_kk(FINAL_INTEGRATE, [](Fix *ifix) {
    ifix->final_integrate();
  });
}

/* ----------------------------------------------------------------------
   fused initial and final integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::fused_integrate(int vflag)
{
  for_each_fix_kk(FINAL_INTEGRATE, [vflag](Fix *ifix) {
    ifix->fused_integrate(vflag);
  });
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyKokkos::end_of_step()
{
  const bigint ntimestep = update->ntimestep;
  for (std::size_t i = 0; i < fixes.size(); ++i)
    if ((callmask[i] & END_OF_STEP) && (ntimestep % fixes[i]->nevery == 0))
      sync_call(fixes[i], [](Fix *ifix) {
        ifix->end_of_step();
      });
}

/* ----------------------------------------------------------------------
   coupling energy call, only for relevant fixes
   each thermostat fix returns this via compute_scalar()
   ecouple = cumulative energy added to reservoir by thermostatting
------------------------------------------------------------------------- */

double ModifyKokkos::energy_couple()
{
  double energy = 0.0;
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->ecouple_flag)
      sync_call(fixes[i], [&energy](Fix *ifix) {
        energy += ifix->compute_scalar();
      });
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
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->energy_global_flag && fixes[i]->thermo_energy)
      sync_call(fixes[i], [&energy](Fix *ifix) {
        energy += ifix->compute_scalar();
      });
  return energy;
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void ModifyKokkos::post_run()
{
  for (std::size_t i = 0; i < fixes.size(); i++) {
    Fix *ifix = fixes[i];
    atomKK->sync(ifix->execution_space, ifix->datamask_read);
    ifix->post_run();
    atomKK->modified(ifix->execution_space, ifix->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force_respa(int vflag, int ilevel)
{
  for_each_fix_kk(PRE_FORCE_RESPA, [=](Fix *ifix) {
    ifix->setup_pre_force_respa(vflag, ilevel);
  });
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  for_each_fix_kk(INITIAL_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->initial_integrate_respa(vflag, ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate_respa(int ilevel, int iloop)
{
  for_each_fix_kk(POST_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->post_integrate_respa(ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force_respa(int vflag, int ilevel, int iloop)
{
  for_each_fix_kk(PRE_FORCE_RESPA, [=](Fix *ifix) {
    ifix->pre_force_respa(vflag, ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
------------------------------------------------------------------------- */

void ModifyKokkos::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (n_post_force_group)
    for_each_fix_kk(GROUP_FIX, [=](Fix *ifix) {
      ifix->post_force_respa(vflag, ilevel, iloop);
    });
  if (n_post_force_respa)
    for_each_fix_kk(POST_FORCE_RESPA, [=](Fix *ifix) {
      ifix->post_force_respa(vflag, ilevel, iloop);
    });
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate_respa(int ilevel, int iloop)
{
  for_each_fix_kk(FINAL_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->final_integrate_respa(ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_exchange()
{
  for_each_fix_kk(MIN_PRE_EXCHANGE, [](Fix *ifix) {
    ifix->min_pre_exchange();
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_neighbor()
{
  for_each_fix_kk(MIN_PRE_NEIGHBOR, [](Fix *ifix) {
    ifix->min_pre_neighbor();
  });
}

/* ----------------------------------------------------------------------
   minimizer post-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_post_neighbor()
{
  for_each_fix_kk(MIN_POST_NEIGHBOR, [](Fix *ifix) {
    ifix->min_post_neighbor();
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_force(int vflag)
{
  for_each_fix_kk(MIN_PRE_FORCE, [vflag](Fix *ifix) {
    ifix->min_pre_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_reverse(int eflag, int vflag)
{
  for_each_fix_kk(MIN_PRE_REVERSE, [=](Fix *ifix) {
    ifix->min_pre_reverse(eflag, vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_post_force(int vflag)
{
  for_each_fix_kk(MIN_POST_FORCE, [vflag](Fix *ifix) {
    ifix->min_post_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double ModifyKokkos::min_energy(double *fextra)
{
  int index = 0;
  double eng = 0.0;
  for_each_fix_kk(MIN_ENERGY, [&](Fix *ifix) {
    eng += ifix->min_energy(&fextra[index]);
    index += ifix->min_dof();
  });
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_store()
{
  for_each_fix_kk(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_store();
  });
}

/* ----------------------------------------------------------------------
   manage state of extra minimizer dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_clearstore()
{
  for_each_fix_kk(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_clearstore();
  });
}

void ModifyKokkos::min_pushstore()
{
  for_each_fix_kk(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_pushstore();
  });
}

void ModifyKokkos::min_popstore()
{
  for_each_fix_kk(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_popstore();
  });
}

/* ----------------------------------------------------------------------
   displace extra minimizer dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_step(double alpha, double *hextra)
{
  int index = 0;
  for_each_fix_kk(MIN_ENERGY, [&](Fix *ifix) {
    ifix->min_step(alpha, &hextra[index]);
    index += ifix->min_dof();
  });
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double ModifyKokkos::max_alpha(double *hextra)
{
  int index = 0;
  double alpha = BIG;
  for_each_fix_kk(MIN_ENERGY, [&](Fix *ifix) {
    double alpha_one = ifix->max_alpha(&hextra[index]);
    alpha = MIN(alpha, alpha_one);
    index += ifix->min_dof();
  });
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_dof()
{
  int ndof = 0;
  for_each_fix_kk(MIN_ENERGY, [&](Fix *ifix) {
    ndof += ifix->min_dof();
  });
  return ndof;
}

/* ----------------------------------------------------------------------
   reset minimizer reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_reset_ref()
{
  int itmpall = 0;
  for_each_fix_kk(MIN_ENERGY, [&](Fix *ifix) {
    if (ifix->min_reset_ref()) itmpall = 1;
  });
  return itmpall;
}
