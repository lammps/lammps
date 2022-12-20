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

#include "modify.h"
#include "style_compute.h"    // IWYU pragma: keep
#include "style_fix.h"        // IWYU pragma: keep

#include "atom.h"
#include "comm.h"
#include "compute.h"    // IWYU pragma: keep
#include "domain.h"
#include "error.h"
#include "fix.h"    // IWYU pragma: keep
#include "group.h"
#include "input.h"
#include "memory.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 4
#define BIG 1.0e20

// template for factory function:
// there will be one instance for each style keyword in the respective style_xxx.h files

template <typename S, typename T> static S *style_creator(LAMMPS *lmp, int narg, char **arg)
{
  return new T(lmp, narg, arg);
}

/* ---------------------------------------------------------------------- */

Modify::Modify(LAMMPS *lmp) : Pointers(lmp)
{
  nfix = maxfix = 0;
  n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = n_post_neighbor = 0;
  n_pre_force = n_pre_reverse = n_post_force_any = 0;
  n_final_integrate = n_end_of_step = 0;
  n_energy_couple = n_energy_global = n_energy_atom = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa_any = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_force = n_min_pre_reverse = 0;
  n_min_post_force = n_min_energy = 0;

  n_timeflag = -1;

  fix = nullptr;
  fmask = nullptr;
  list_initial_integrate = list_post_integrate = nullptr;
  list_pre_exchange = list_pre_neighbor = list_post_neighbor = nullptr;
  list_pre_force = list_pre_reverse = nullptr;
  list_post_force = list_post_force_group = nullptr;
  list_final_integrate = list_end_of_step = nullptr;
  list_energy_couple = list_energy_global = list_energy_atom = nullptr;
  list_initial_integrate_respa = list_post_integrate_respa = nullptr;
  list_pre_force_respa = list_post_force_respa = nullptr;
  list_final_integrate_respa = nullptr;
  list_min_pre_exchange = list_min_pre_neighbor = list_min_post_neighbor = nullptr;
  list_min_pre_force = list_min_pre_reverse = list_min_post_force = nullptr;
  list_min_energy = nullptr;

  end_of_step_every = nullptr;

  list_timeflag = nullptr;

  nfix_restart_global = 0;
  id_restart_global = style_restart_global = nullptr;
  state_restart_global = nullptr;
  used_restart_global = nullptr;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = nullptr;
  index_restart_peratom = used_restart_peratom = nullptr;

  ncompute = maxcompute = 0;
  compute = nullptr;

  create_factories();
}

void _noopt Modify::create_factories()
{
  // fill map with fixes listed in style_fix.h

  fix_map = new FixCreatorMap();

#define FIX_CLASS
#define FixStyle(key, Class) (*fix_map)[#key] = &style_creator<Fix, Class>;
#include "style_fix.h"    // IWYU pragma: keep
#undef FixStyle
#undef FIX_CLASS

  // fill map with computes listed in style_compute.h

  compute_map = new ComputeCreatorMap();

#define COMPUTE_CLASS
#define ComputeStyle(key, Class) (*compute_map)[#key] = &style_creator<Compute, Class>;
#include "style_compute.h"    // IWYU pragma: keep
#undef ComputeStyle
#undef COMPUTE_CLASS
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes
  // do it via delete_fix() so callbacks in Atom are also updated correctly

  while (nfix) delete_fix(0);
  memory->sfree(fix);
  memory->destroy(fmask);

  // delete all computes

  for (int i = 0; i < ncompute; i++) delete compute[i];
  memory->sfree(compute);

  delete[] list_initial_integrate;
  delete[] list_post_integrate;
  delete[] list_pre_exchange;
  delete[] list_pre_neighbor;
  delete[] list_post_neighbor;
  delete[] list_pre_force;
  delete[] list_pre_reverse;
  delete[] list_post_force;
  delete[] list_post_force_group;
  delete[] list_final_integrate;
  delete[] list_end_of_step;
  delete[] list_energy_couple;
  delete[] list_energy_global;
  delete[] list_energy_atom;
  delete[] list_initial_integrate_respa;
  delete[] list_post_integrate_respa;
  delete[] list_pre_force_respa;
  delete[] list_post_force_respa;
  delete[] list_final_integrate_respa;
  delete[] list_min_pre_exchange;
  delete[] list_min_pre_neighbor;
  delete[] list_min_post_neighbor;
  delete[] list_min_pre_force;
  delete[] list_min_pre_reverse;
  delete[] list_min_post_force;
  delete[] list_min_energy;

  delete[] end_of_step_every;
  delete[] list_timeflag;

  restart_deallocate(0);

  delete compute_map;
  delete fix_map;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i, j;

  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate(1);

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes

  for (i = 0; i < ncompute; i++) {
    compute[i]->init();
    compute[i]->invoked_scalar = -1;
    compute[i]->invoked_vector = -1;
    compute[i]->invoked_array = -1;
    compute[i]->invoked_peratom = -1;
    compute[i]->invoked_local = -1;
  }
  addstep_compute_all(update->ntimestep);

  // init each fix
  // should not need to come before compute init
  //   used to b/c temperature computes called fix->dof() in their init,
  //   and fix rigid required its own init before its dof() could be called,
  //   but computes now do their DOF in setup()

  for (i = 0; i < nfix; i++) fix[i]->init();

  // set global flag if any fix has its restart_pbc flag set

  restart_pbc_any = 0;
  for (i = 0; i < nfix; i++)
    if (fix[i]->restart_pbc) restart_pbc_any = 1;

  // create lists of fixes to call at each stage of run
  // needs to happen after init() of computes
  //   b/c a compute::init() can delete a fix, e.g. compute chunk/atom

  list_init(INITIAL_INTEGRATE, n_initial_integrate, list_initial_integrate);
  list_init(POST_INTEGRATE, n_post_integrate, list_post_integrate);
  list_init(PRE_EXCHANGE, n_pre_exchange, list_pre_exchange);
  list_init(PRE_NEIGHBOR, n_pre_neighbor, list_pre_neighbor);
  list_init(POST_NEIGHBOR, n_post_neighbor, list_post_neighbor);
  list_init(PRE_FORCE, n_pre_force, list_pre_force);
  list_init(PRE_REVERSE, n_pre_reverse, list_pre_reverse);
  list_init(POST_FORCE, n_post_force, list_post_force);
  list_init_post_force_group(n_post_force_group, list_post_force_group);
  list_init(FINAL_INTEGRATE, n_final_integrate, list_final_integrate);
  list_init_end_of_step(END_OF_STEP, n_end_of_step, list_end_of_step);
  list_init_energy_couple(n_energy_couple, list_energy_couple);
  list_init_energy_global(n_energy_global, list_energy_global);
  list_init_energy_atom(n_energy_atom, list_energy_atom);

  list_init(INITIAL_INTEGRATE_RESPA, n_initial_integrate_respa, list_initial_integrate_respa);
  list_init(POST_INTEGRATE_RESPA, n_post_integrate_respa, list_post_integrate_respa);
  list_init(POST_FORCE_RESPA, n_post_force_respa, list_post_force_respa);
  list_init(PRE_FORCE_RESPA, n_pre_force_respa, list_pre_force_respa);
  list_init(FINAL_INTEGRATE_RESPA, n_final_integrate_respa, list_final_integrate_respa);

  list_init(MIN_PRE_EXCHANGE, n_min_pre_exchange, list_min_pre_exchange);
  list_init(MIN_PRE_NEIGHBOR, n_min_pre_neighbor, list_min_pre_neighbor);
  list_init(MIN_POST_NEIGHBOR, n_min_post_neighbor, list_min_post_neighbor);
  list_init(MIN_PRE_FORCE, n_min_pre_force, list_min_pre_force);
  list_init(MIN_PRE_REVERSE, n_min_pre_reverse, list_min_pre_reverse);
  list_init(MIN_POST_FORCE, n_min_post_force, list_min_post_force);
  list_init(MIN_ENERGY, n_min_energy, list_min_energy);

  // two post_force_any counters used by integrators add in post_force_group

  n_post_force_any = n_post_force + n_post_force_group;
  n_post_force_respa_any = n_post_force_respa + n_post_force_group;

  // create list of computes that store invocation times

  list_init_compute();

  // error if any fix or compute is using a dynamic group when not allowed

  for (i = 0; i < nfix; i++)
    if (!fix[i]->dynamic_group_allow && group->dynamic[fix[i]->igroup])
      error->all(FLERR, "Fix {} does not allow use with a dynamic group", fix[i]->style);

  for (i = 0; i < ncompute; i++)
    if (!compute[i]->dynamic_group_allow && group->dynamic[compute[i]->igroup])
      error->all(FLERR, "Compute {} does not allow use with a dynamic group", compute[i]->style);

  // warn if any particle is time integrated more than once

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++) flag[i] = 0;

  int groupbit;
  for (i = 0; i < nfix; i++) {
    if (fix[i]->time_integrate == 0) continue;
    groupbit = fix[i]->groupbit;
    for (j = 0; j < nlocal; j++)
      if (mask[j] & groupbit) flag[j]++;
  }

  int check = 0;
  for (i = 0; i < nlocal; i++)
    if (flag[i] > 1) check = 1;

  delete[] flag;

  int checkall;
  MPI_Allreduce(&check, &checkall, 1, MPI_INT, MPI_SUM, world);
  if (comm->me == 0 && checkall)
    error->warning(FLERR, "One or more atoms are time integrated more than once");
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup(int vflag)
{
  // compute setup needs to come before fix setup
  //   b/c NH fixes need DOF of temperature computes
  // fix group setup() is special case since populates a dynamic group
  //   needs to be done before temperature compute setup

  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0) fix[i]->setup(vflag);

  for (int i = 0; i < ncompute; i++) compute[i]->setup();

  if (update->whichflag == 1)
    for (int i = 0; i < nfix; i++) fix[i]->setup(vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < nfix; i++) fix[i]->min_setup(vflag);
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void Modify::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    for (int i = 0; i < n_pre_exchange; i++) fix[list_pre_exchange[i]]->setup_pre_exchange();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_exchange; i++)
      fix[list_min_pre_exchange[i]]->setup_pre_exchange();
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_neighbor; i++) fix[list_pre_neighbor[i]]->setup_pre_neighbor();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_neighbor; i++)
      fix[list_min_pre_neighbor[i]]->setup_pre_neighbor();
}

/* ----------------------------------------------------------------------
   setup post_neighbor call, only for fixes that define post_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_post_neighbor()
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_post_neighbor; i++) fix[list_post_neighbor[i]]->setup_post_neighbor();
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_post_neighbor; i++)
      fix[list_min_post_neighbor[i]]->setup_post_neighbor();
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_force(int vflag)
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_force; i++) fix[list_pre_force[i]]->setup_pre_force(vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_force; i++) fix[list_min_pre_force[i]]->setup_pre_force(vflag);
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_reverse(int eflag, int vflag)
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_reverse; i++)
      fix[list_pre_reverse[i]]->setup_pre_reverse(eflag, vflag);
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_reverse; i++)
      fix[list_min_pre_reverse[i]]->setup_pre_reverse(eflag, vflag);
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate(int vflag)
{
  for (int i = 0; i < n_initial_integrate; i++)
    fix[list_initial_integrate[i]]->initial_integrate(vflag);
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate()
{
  for (int i = 0; i < n_post_integrate; i++) fix[list_post_integrate[i]]->post_integrate();
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_exchange()
{
  for (int i = 0; i < n_pre_exchange; i++) fix[list_pre_exchange[i]]->pre_exchange();
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_neighbor()
{
  for (int i = 0; i < n_pre_neighbor; i++) fix[list_pre_neighbor[i]]->pre_neighbor();
}

/* ----------------------------------------------------------------------
   post_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_neighbor()
{
  for (int i = 0; i < n_post_neighbor; i++) fix[list_post_neighbor[i]]->post_neighbor();
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force(int vflag)
{
  for (int i = 0; i < n_pre_force; i++) fix[list_pre_force[i]]->pre_force(vflag);
}
/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_reverse(int eflag, int vflag)
{
  for (int i = 0; i < n_pre_reverse; i++) fix[list_pre_reverse[i]]->pre_reverse(eflag, vflag);
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
     they are not in n_post_force count
------------------------------------------------------------------------- */

void Modify::post_force(int vflag)
{
  if (n_post_force_group) {
    for (int i = 0; i < n_post_force_group; i++) fix[list_post_force_group[i]]->post_force(vflag);
  }

  if (n_post_force) {
    for (int i = 0; i < n_post_force; i++) fix[list_post_force[i]]->post_force(vflag);
  }
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++) fix[list_final_integrate[i]]->final_integrate();
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0) fix[list_end_of_step[i]]->end_of_step();
}

/* ----------------------------------------------------------------------
   coupling energy call, only for relevant fixes
   each thermostsat fix returns this via compute_scalar()
   ecouple = cumulative energy added to reservoir by thermostatting
------------------------------------------------------------------------- */

double Modify::energy_couple()
{
  double energy = 0.0;
  for (int i = 0; i < n_energy_couple; i++) energy += fix[list_energy_couple[i]]->compute_scalar();
  return energy;
}

/* ----------------------------------------------------------------------
   global energy call, only for relevant fixes
   they return energy via compute_scalar()
   called by compute pe
------------------------------------------------------------------------- */

double Modify::energy_global()
{
  double energy = 0.0;
  for (int i = 0; i < n_energy_global; i++) energy += fix[list_energy_global[i]]->compute_scalar();
  return energy;
}

/* ----------------------------------------------------------------------
   peratom energy call, only for relevant fixes
   called by compute pe/atom
------------------------------------------------------------------------- */

void Modify::energy_atom(int nlocal, double *energy)
{
  int i, j;
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

void Modify::post_run()
{
  for (int i = 0; i < nfix; i++) fix[i]->post_run();

  // must reset this to its default value, since computes may be added
  // or removed between runs and with this change we will redirect any
  // calls to addstep_compute() to addstep_compute_all() instead.
  n_timeflag = -1;
}

/* ----------------------------------------------------------------------
   create_attribute call
   invoked when an atom is added to system during a run
   necessary so that fixes and computes that store per-atom
     state can initialize that state for the new atom N
   computes can store per-atom state via a fix like fix STORE
     compute has the create_attribute flag, not fix STORE
------------------------------------------------------------------------- */

void Modify::create_attribute(int n)
{
  for (int i = 0; i < nfix; i++)
    if (fix[i]->create_attribute) fix[i]->set_arrays(n);
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->create_attribute) compute[i]->set_arrays(n);
  input->variable->set_arrays(n);
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::setup_pre_force_respa(int vflag, int ilevel)
{
  for (int i = 0; i < n_pre_force_respa; i++)
    fix[list_pre_force_respa[i]]->setup_pre_force_respa(vflag, ilevel);
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_initial_integrate_respa; i++)
    fix[list_initial_integrate_respa[i]]->initial_integrate_respa(vflag, ilevel, iloop);
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_post_integrate_respa; i++)
    fix[list_post_integrate_respa[i]]->post_integrate_respa(ilevel, iloop);
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_pre_force_respa; i++)
    fix[list_pre_force_respa[i]]->pre_force_respa(vflag, ilevel, iloop);
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
------------------------------------------------------------------------- */

void Modify::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (n_post_force_group) {
    for (int i = 0; i < n_post_force_group; i++)
      fix[list_post_force_group[i]]->post_force_respa(vflag, ilevel, iloop);
  }

  if (n_post_force_respa) {
    for (int i = 0; i < n_post_force_respa; i++)
      fix[list_post_force_respa[i]]->post_force_respa(vflag, ilevel, iloop);
  }
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_final_integrate_respa; i++)
    fix[list_final_integrate_respa[i]]->final_integrate_respa(ilevel, iloop);
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_exchange()
{
  for (int i = 0; i < n_min_pre_exchange; i++) fix[list_min_pre_exchange[i]]->min_pre_exchange();
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_neighbor()
{
  for (int i = 0; i < n_min_pre_neighbor; i++) fix[list_min_pre_neighbor[i]]->min_pre_neighbor();
}

/* ----------------------------------------------------------------------
   minimizer post-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_neighbor()
{
  for (int i = 0; i < n_min_post_neighbor; i++) fix[list_min_post_neighbor[i]]->min_post_neighbor();
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_force(int vflag)
{
  for (int i = 0; i < n_min_pre_force; i++) fix[list_min_pre_force[i]]->min_pre_force(vflag);
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_reverse(int eflag, int vflag)
{
  for (int i = 0; i < n_min_pre_reverse; i++)
    fix[list_min_pre_reverse[i]]->min_pre_reverse(eflag, vflag);
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_force(int vflag)
{
  for (int i = 0; i < n_min_post_force; i++) fix[list_min_post_force[i]]->min_post_force(vflag);
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double Modify::min_energy(double *fextra)
{
  int ifix, index;

  index = 0;
  double eng = 0.0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    eng += fix[ifix]->min_energy(&fextra[index]);
    index += fix[ifix]->min_dof();
  }
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_store()
{
  for (int i = 0; i < n_min_energy; i++) fix[list_min_energy[i]]->min_store();
}

/* ----------------------------------------------------------------------
   manage state of extra minimizer dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_clearstore()
{
  for (int i = 0; i < n_min_energy; i++) fix[list_min_energy[i]]->min_clearstore();
}

void Modify::min_pushstore()
{
  for (int i = 0; i < n_min_energy; i++) fix[list_min_energy[i]]->min_pushstore();
}

void Modify::min_popstore()
{
  for (int i = 0; i < n_min_energy; i++) fix[list_min_energy[i]]->min_popstore();
}

/* ----------------------------------------------------------------------
   displace extra minimizer dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_step(double alpha, double *hextra)
{
  int ifix, index;

  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    fix[ifix]->min_step(alpha, &hextra[index]);
    index += fix[ifix]->min_dof();
  }
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double Modify::max_alpha(double *hextra)
{
  int ifix, index;

  double alpha = BIG;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    double alpha_one = fix[ifix]->max_alpha(&hextra[index]);
    alpha = MIN(alpha, alpha_one);
    index += fix[ifix]->min_dof();
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_dof()
{
  int ndof = 0;
  for (int i = 0; i < n_min_energy; i++) ndof += fix[list_min_energy[i]]->min_dof();
  return ndof;
}

/* ----------------------------------------------------------------------
   reset minimizer reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_reset_ref()
{
  int itmp, itmpall;
  itmpall = 0;
  for (int i = 0; i < n_min_energy; i++) {
    itmp = fix[list_min_energy[i]]->min_reset_ref();
    if (itmp) itmpall = 1;
  }
  return itmpall;
}

/* ----------------------------------------------------------------------
   reset grids for any Fix or Compute that uses distributed grids
   called by load balancer when proc sub-domains change
------------------------------------------------------------------------- */

void Modify::reset_grid()
{
  for (int i = 0; i < nfix; i++)
    if (fix[i]->pergrid_flag) fix[i]->reset_grid();
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->pergrid_flag) compute[i]->reset_grid();
}

/* ----------------------------------------------------------------------
   add a new fix or replace one with same ID
------------------------------------------------------------------------- */

Fix *Modify::add_fix(int narg, char **arg, int trysuffix)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "fix", error);

  // cannot define fix before box exists unless style is in exception list
  // don't like this way of checking for exceptions by adding fixes to list,
  //   but can't think of better way
  // too late if instantiate fix, then check flag set in fix constructor,
  //   since some fixes access domain settings in their constructor
  // nullptr must be last entry in this list

  // clang-format off
  const char *exceptions[] =
    {"GPU", "OMP", "INTEL", "property/atom", "cmap", "cmap3", "rx",
     "deprecated", "STORE/KIM", "amoeba/pitorsion", "amoeba/bitorsion",
     nullptr};
  // clang-format on

  if (domain->box_exist == 0) {
    int m;
    for (m = 0; exceptions[m] != nullptr; m++)
      if (strcmp(arg[2], exceptions[m]) == 0) break;
    if (exceptions[m] == nullptr) error->all(FLERR, "Fix command before simulation box is defined");
  }

  // check group ID

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR, "Could not find fix group ID {}", arg[1]);

  // if fix ID exists:
  //   set newflag = 0 so create new fix in same location in fix list
  //   error if new style does not match old style
  //     since can't replace it (all when-to-invoke ptrs would be invalid)
  //   warn if new group != old group
  //   delete old fix, but do not call update_callback(),
  //     since will replace this fix and thus other fix locs will not change
  //   set ptr to a null pointer in case new fix scans list of fixes,
  //     e.g. scan will occur in add_callback() if called by new fix
  // if fix ID does not exist:
  //   set newflag = 1 so create new fix
  //   extend fix and fmask lists as necessary

  int ifix, newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0], fix[ifix]->id) == 0) break;

  if (ifix < nfix) {
    newflag = 0;

    int match = 0;
    if (strcmp(arg[2], fix[ifix]->style) == 0) match = 1;
    if (!match && trysuffix && lmp->suffix_enable) {
      if (lmp->non_pair_suffix()) {
        std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
        if (estyle == fix[ifix]->style) match = 1;
      }
      if (lmp->suffix2) {
        std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
        if (estyle == fix[ifix]->style) match = 1;
      }
    }
    if (!match) error->all(FLERR, "Replacing a fix, but new style != old style");

    if (fix[ifix]->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Replacing a fix, but new group != old group");
    delete fix[ifix];
    fix[ifix] = nullptr;

  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix, maxfix * sizeof(Fix *), "modify:fix");
      memory->grow(fmask, maxfix, "modify:fmask");
    }
  }

  // create the Fix
  // try first with suffix appended

  fix[ifix] = nullptr;

  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
      if (fix_map->find(estyle) != fix_map->end()) {
        FixCreator &fix_creator = (*fix_map)[estyle];
        fix[ifix] = fix_creator(lmp, narg, arg);
        delete[] fix[ifix]->style;
        fix[ifix]->style = utils::strdup(estyle);
      }
    }
    if ((fix[ifix] == nullptr) && lmp->suffix2) {
      std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
      if (fix_map->find(estyle) != fix_map->end()) {
        FixCreator &fix_creator = (*fix_map)[estyle];
        fix[ifix] = fix_creator(lmp, narg, arg);
        delete[] fix[ifix]->style;
        fix[ifix]->style = utils::strdup(estyle);
      }
    }
  }

  if ((fix[ifix] == nullptr) && (fix_map->find(arg[2]) != fix_map->end())) {
    FixCreator &fix_creator = (*fix_map)[arg[2]];
    fix[ifix] = fix_creator(lmp, narg, arg);
  }

  if (fix[ifix] == nullptr) error->all(FLERR, utils::check_packages_for_style("fix", arg[2], lmp));

  // increment nfix and update fix_list vector (if new)

  if (newflag) {
    nfix++;
    fix_list = std::vector<Fix *>(fix, fix + nfix);
  }

  // post_constructor() can call virtual methods in parent or child
  //   which would otherwise not yet be visible in child class
  // post_constructor() allows new fix to create other fixes
  // nfix increment must come first so recursive call to add_fix within
  //   post_constructor() will see updated nfix

  fix[ifix]->post_constructor();

  // check if Fix is in restart_global list
  // if yes, pass state info to the Fix so it can reset itself

  for (int i = 0; i < nfix_restart_global; i++)
    if ((strcmp(id_restart_global[i], fix[ifix]->id) == 0) &&
        (utils::strip_style_suffix(fix[ifix]->style, lmp) == style_restart_global[i])) {
      fix[ifix]->restart(state_restart_global[i]);
      used_restart_global[i] = 1;
      fix[ifix]->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting global fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       fix[ifix]->style, fix[ifix]->id);
    }

  // check if Fix is in restart_peratom list
  // if yes, loop over atoms so they can extract info from atom->extra array

  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strcmp(id_restart_peratom[i], fix[ifix]->id) == 0 &&
        strcmp(style_restart_peratom[i], fix[ifix]->style) == 0) {
      used_restart_peratom[i] = 1;
      for (int j = 0; j < atom->nlocal; j++) fix[ifix]->unpack_restart(j, index_restart_peratom[i]);
      fix[ifix]->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting peratom fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       fix[ifix]->style, fix[ifix]->id);
    }

  // set fix mask values

  fmask[ifix] = fix[ifix]->setmask();

  // return pointer to fix

  return fix[ifix];
}

/* ----------------------------------------------------------------------
   convenience function to allow adding a fix from a single string
------------------------------------------------------------------------- */

Fix *Modify::add_fix(const std::string &fixcmd, int trysuffix)
{
  auto args = utils::split_words(fixcmd);
  std::vector<char *> newarg(args.size());
  int i = 0;
  for (const auto &arg : args) { newarg[i++] = (char *) arg.c_str(); }
  return add_fix(args.size(), newarg.data(), trysuffix);
}

/* ----------------------------------------------------------------------
   replace replaceID fix with a new fix
   this is used by callers to preserve ordering of fixes
   e.g. create replaceID as a FixDummy instance early in the input script
        replace it later with the desired Fix instance
------------------------------------------------------------------------- */

Fix *Modify::replace_fix(const char *replaceID, int narg, char **arg, int trysuffix)
{
  auto oldfix = get_fix_by_id(replaceID);
  if (!oldfix) error->all(FLERR, "Modify replace_fix ID {} could not be found", replaceID);

  // change ID, igroup, style of fix being replaced to match new fix
  // requires some error checking on arguments for new fix

  if (narg < 3) error->all(FLERR, "Not enough arguments for replace_fix invocation");
  if (get_fix_by_id(arg[0])) error->all(FLERR, "Replace_fix ID {} is already in use", arg[0]);

  delete[] oldfix->id;
  oldfix->id = utils::strdup(arg[0]);

  int jgroup = group->find(arg[1]);
  if (jgroup == -1) error->all(FLERR, "Could not find replace_fix group ID {}", arg[1]);
  oldfix->igroup = jgroup;

  delete[] oldfix->style;
  oldfix->style = utils::strdup(arg[2]);

  // invoke add_fix
  // it will find and overwrite the replaceID fix

  return add_fix(narg, arg, trysuffix);
}

/* ----------------------------------------------------------------------
   convenience function to allow replacing a fix from a single string
------------------------------------------------------------------------- */

Fix *Modify::replace_fix(const std::string &oldfix, const std::string &fixcmd, int trysuffix)
{
  auto args = utils::split_words(fixcmd);
  std::vector<char *> newarg(args.size());
  int i = 0;
  for (const auto &arg : args) { newarg[i++] = (char *) arg.c_str(); }
  return replace_fix(oldfix.c_str(), args.size(), newarg.data(), trysuffix);
}

/* ----------------------------------------------------------------------
   modify a Fix's parameters
------------------------------------------------------------------------- */

void Modify::modify_fix(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);

  auto ifix = get_fix_by_id(arg[0]);
  if (!ifix) error->all(FLERR, "Could not find fix_modify ID {}", arg[0]);
  ifix->modify_params(narg - 1, &arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
   Atom class must update indices in its list of callbacks to fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(const std::string &id)
{
  int ifix = find_fix(id);
  if (ifix < 0) error->all(FLERR, "Could not find fix ID {} to delete", id);
  delete_fix(ifix);
}

void Modify::delete_fix(int ifix)
{
  if ((ifix < 0) || (ifix >= nfix)) return;

  // delete instance and move other Fixes and fmask down in list one slot

  delete fix[ifix];
  atom->update_callback(ifix);

  for (int i = ifix + 1; i < nfix; i++) fix[i - 1] = fix[i];
  for (int i = ifix + 1; i < nfix; i++) fmask[i - 1] = fmask[i];
  nfix--;
  fix_list = std::vector<Fix *>(fix, fix + nfix);
}

/* ----------------------------------------------------------------------
   find a fix by ID
   return index of fix or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_fix(const std::string &id)
{
  if (id.empty()) return -1;
  for (int ifix = 0; ifix < nfix; ifix++)
    if (id == fix[ifix]->id) return ifix;
  return -1;
}

/* ----------------------------------------------------------------------
   look up pointer to Fix class by fix-ID
   return null pointer if ID not found
------------------------------------------------------------------------- */

Fix *Modify::get_fix_by_id(const std::string &id) const
{
  if (id.empty()) return nullptr;
  for (int ifix = 0; ifix < nfix; ifix++)
    if (id == fix[ifix]->id) return fix[ifix];
  return nullptr;
}

/* ----------------------------------------------------------------------
   look up pointer to fixes by fix style name
   return vector of matching pointers
------------------------------------------------------------------------- */

const std::vector<Fix *> Modify::get_fix_by_style(const std::string &style) const
{
  std::vector<Fix *> matches;
  if (style.empty()) return matches;

  for (int ifix = 0; ifix < nfix; ifix++)
    if (utils::strmatch(fix[ifix]->style, style)) matches.push_back(fix[ifix]);

  return matches;
}

/* ----------------------------------------------------------------------
   return list of fixes as vector
------------------------------------------------------------------------- */

const std::vector<Fix *> &Modify::get_fix_list()
{
  fix_list = std::vector<Fix *>(fix, fix + nfix);
  return fix_list;
}

/* ----------------------------------------------------------------------
   check for fix associated with package name in compiled list
   return 1 if found else 0
   used to determine whether LAMMPS was built with
     GPU, INTEL, OPENMP packages, which have their own fixes
------------------------------------------------------------------------- */

int Modify::check_package(const char *package_fix_name)
{
  if (fix_map->find(package_fix_name) == fix_map->end()) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check if the group indicated by groupbit overlaps with any
   currently existing rigid fixes. return 1 in this case otherwise 0
------------------------------------------------------------------------- */

int Modify::check_rigid_group_overlap(int groupbit)
{
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;
  int dim;

  int n = 0;
  for (int ifix = 0; ifix < nfix; ifix++) {
    if (utils::strmatch(fix[ifix]->style, "^rigid")) {
      const int *const body = (const int *) fix[ifix]->extract("body", dim);
      if ((body == nullptr) || (dim != 1)) break;

      for (int i = 0; (i < nlocal) && (n == 0); ++i)
        if ((mask[i] & groupbit) && (body[i] >= 0)) ++n;
    }
  }

  int n_all = 0;
  MPI_Allreduce(&n, &n_all, 1, MPI_INT, MPI_SUM, world);

  if (n_all > 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if the atoms in the group indicated by groupbit _and_ region
   indicated by regionid overlap with any currently existing rigid fixes.
   return 1 in this case, otherwise 0
------------------------------------------------------------------------- */

int Modify::check_rigid_region_overlap(int groupbit, Region *reg)
{
  const int *const mask = atom->mask;
  const double *const *const x = atom->x;
  const int nlocal = atom->nlocal;
  int dim;

  int n = 0;
  reg->prematch();
  for (int ifix = 0; ifix < nfix; ifix++) {
    if (strncmp("rigid", fix[ifix]->style, 5) == 0) {
      const int *const body = (const int *) fix[ifix]->extract("body", dim);
      if ((body == nullptr) || (dim != 1)) break;

      for (int i = 0; (i < nlocal) && (n == 0); ++i)
        if ((mask[i] & groupbit) && (body[i] >= 0) && reg->match(x[i][0], x[i][1], x[i][2])) ++n;
    }
  }

  int n_all = 0;
  MPI_Allreduce(&n, &n_all, 1, MPI_INT, MPI_SUM, world);

  if (n_all > 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if the atoms in the selection list (length atom->nlocal,
   content: 1 if atom is contained, 0 if not) overlap with currently
   existing rigid fixes. return 1 in this case otherwise 0
------------------------------------------------------------------------- */

int Modify::check_rigid_list_overlap(int *select)
{
  const int nlocal = atom->nlocal;
  int dim;

  int n = 0;
  for (int ifix = 0; ifix < nfix; ifix++) {
    if (utils::strmatch(fix[ifix]->style, "^rigid")) {
      const int *const body = (const int *) fix[ifix]->extract("body", dim);
      if ((body == nullptr) || (dim != 1)) break;

      for (int i = 0; (i < nlocal) && (n == 0); ++i)
        if ((body[i] >= 0) && select[i]) ++n;
    }
  }

  int n_all = 0;
  MPI_Allreduce(&n, &n_all, 1, MPI_INT, MPI_SUM, world);

  if (n_all > 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   add a new compute
------------------------------------------------------------------------- */

Compute *Modify::add_compute(int narg, char **arg, int trysuffix)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "compute", error);

  // error check

  if (get_compute_by_id(arg[0])) error->all(FLERR, "Reuse of compute ID '{}'", arg[0]);

  // extend Compute list if necessary

  if (ncompute == maxcompute) {
    maxcompute += DELTA;
    compute =
        (Compute **) memory->srealloc(compute, maxcompute * sizeof(Compute *), "modify:compute");
  }

  // create the Compute
  // try first with suffix appended

  compute[ncompute] = nullptr;

  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
      if (compute_map->find(estyle) != compute_map->end()) {
        ComputeCreator &compute_creator = (*compute_map)[estyle];
        compute[ncompute] = compute_creator(lmp, narg, arg);
        delete[] compute[ncompute]->style;
        compute[ncompute]->style = utils::strdup(estyle);
      }
    }
    if (compute[ncompute] == nullptr && lmp->suffix2) {
      std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
      if (compute_map->find(estyle) != compute_map->end()) {
        ComputeCreator &compute_creator = (*compute_map)[estyle];
        compute[ncompute] = compute_creator(lmp, narg, arg);
        delete[] compute[ncompute]->style;
        compute[ncompute]->style = utils::strdup(estyle);
      }
    }
  }

  if (compute[ncompute] == nullptr && compute_map->find(arg[2]) != compute_map->end()) {
    ComputeCreator &compute_creator = (*compute_map)[arg[2]];
    compute[ncompute] = compute_creator(lmp, narg, arg);
  }

  if (compute[ncompute] == nullptr)
    error->all(FLERR, utils::check_packages_for_style("compute", arg[2], lmp));

  compute_list = std::vector<Compute *>(compute, compute + ncompute + 1);
  return compute[ncompute++];
}

/* ----------------------------------------------------------------------
   convenience function to allow adding a compute from a single string
------------------------------------------------------------------------- */

Compute *Modify::add_compute(const std::string &computecmd, int trysuffix)
{
  auto args = utils::split_words(computecmd);
  std::vector<char *> newarg(args.size());
  int i = 0;
  for (const auto &arg : args) { newarg[i++] = (char *) arg.c_str(); }
  return add_compute(args.size(), newarg.data(), trysuffix);
}

/* ----------------------------------------------------------------------
   modify a Compute's parameters
------------------------------------------------------------------------- */

void Modify::modify_compute(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "compute_modify", error);

  // lookup Compute ID

  auto icompute = get_compute_by_id(arg[0]);
  if (!icompute) error->all(FLERR, "Could not find compute_modify ID {}", arg[0]);
  icompute->modify_params(narg - 1, &arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Compute from list of Computes
------------------------------------------------------------------------- */

void Modify::delete_compute(const std::string &id)
{
  int icompute = find_compute(id);
  if (icompute < 0) error->all(FLERR, "Could not find compute ID {} to delete", id);
  delete_compute(icompute);
}

void Modify::delete_compute(int icompute)
{
  if ((icompute < 0) || (icompute >= ncompute)) return;

  // delete and move other Computes down in list one slot

  delete compute[icompute];
  for (int i = icompute + 1; i < ncompute; i++) compute[i - 1] = compute[i];
  ncompute--;
  compute_list = std::vector<Compute *>(compute, compute + ncompute);
}

/* ----------------------------------------------------------------------
   find a compute by ID
   return index of compute or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_compute(const std::string &id)
{
  if (id.empty()) return -1;
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (id == compute[icompute]->id) return icompute;
  return -1;
}

/* ----------------------------------------------------------------------
   look up pointer to Compute class by compute-ID
   return null pointer if ID not found
------------------------------------------------------------------------- */

Compute *Modify::get_compute_by_id(const std::string &id) const
{
  if (id.empty()) return nullptr;
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (id == compute[icompute]->id) return compute[icompute];
  return nullptr;
}

/* ----------------------------------------------------------------------
   look up pointers to computes by compute style name
   return vector with matching pointers
------------------------------------------------------------------------- */

const std::vector<Compute *> Modify::get_compute_by_style(const std::string &style) const
{
  std::vector<Compute *> matches;
  if (style.empty()) return matches;

  for (int icompute = 0; icompute < ncompute; icompute++)
    if (utils::strmatch(compute[icompute]->style, style)) matches.push_back(compute[icompute]);

  return matches;
}

/* ----------------------------------------------------------------------
   return vector with Computes
------------------------------------------------------------------------- */

const std::vector<Compute *> &Modify::get_compute_list()
{
  compute_list = std::vector<Compute *>(compute, compute + ncompute);
  return compute_list;
}

/* ----------------------------------------------------------------------
   clear invoked flag of all computes
   called everywhere that computes are used, before computes are invoked
   invoked flag used to avoid re-invoking same compute multiple times
   and to flag computes that store invocation times as having been invoked
------------------------------------------------------------------------- */

void Modify::clearstep_compute()
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    compute[icompute]->invoked_flag = Compute::INVOKED_NONE;
}

/* ----------------------------------------------------------------------
   loop over computes that store invocation times
   if its invoked flag set on this timestep, schedule next invocation
   called everywhere that computes are used, after computes are invoked
------------------------------------------------------------------------- */

void Modify::addstep_compute(bigint newstep)
{
  // If we are called before the first run init, n_timeflag is not yet
  // initialized, thus defer to addstep_compute_all() instead

  if (n_timeflag < 0) {
    addstep_compute_all(newstep);
    return;
  }

  for (int icompute = 0; icompute < n_timeflag; icompute++)
    if (compute[list_timeflag[icompute]]->invoked_flag)
      compute[list_timeflag[icompute]]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   loop over all computes
   schedule next invocation for those that store invocation times
   called when not sure what computes will be needed on newstep
   do not loop only over n_timeflag, since may not be set yet
------------------------------------------------------------------------- */

void Modify::addstep_compute_all(bigint newstep)
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (compute[icompute]->timeflag) compute[icompute]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   write to restart file for all Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
------------------------------------------------------------------------- */

void Modify::write_restart(FILE *fp)
{
  int me = comm->me;

  int count = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_global) count++;

  if (me == 0) fwrite(&count, sizeof(int), 1, fp);

  int n;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_global) {
      if (me == 0) {
        n = strlen(fix[i]->id) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(fix[i]->id, sizeof(char), n, fp);
        auto fix_style = utils::strip_style_suffix(fix[i]->style, lmp);
        n = fix_style.size() + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(fix_style.c_str(), sizeof(char), n, fp);
      }
      fix[i]->write_restart(fp);
    }

  count = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_peratom) count++;

  if (me == 0) fwrite(&count, sizeof(int), 1, fp);

  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_peratom) {
      int maxsize_restart = fix[i]->maxsize_restart();
      if (me == 0) {
        n = strlen(fix[i]->id) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(fix[i]->id, sizeof(char), n, fp);
        n = strlen(fix[i]->style) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(fix[i]->style, sizeof(char), n, fp);
        fwrite(&maxsize_restart, sizeof(int), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   read in restart file data on all previously defined Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
   return maxsize of extra info that will be stored with any atom
------------------------------------------------------------------------- */

int Modify::read_restart(FILE *fp)
{
  // nfix_restart_global = # of restart entries with global state info

  int me = comm->me;
  if (me == 0) utils::sfread(FLERR, &nfix_restart_global, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&nfix_restart_global, 1, MPI_INT, 0, world);

  // allocate space for each entry

  if (nfix_restart_global) {
    id_restart_global = new char *[nfix_restart_global];
    style_restart_global = new char *[nfix_restart_global];
    state_restart_global = new char *[nfix_restart_global];
    used_restart_global = new int[nfix_restart_global];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, chunk of state data

  int n;
  for (int i = 0; i < nfix_restart_global; i++) {
    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    id_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, id_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(id_restart_global[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    style_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, style_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(style_restart_global[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    state_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, state_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(state_restart_global[i], n, MPI_CHAR, 0, world);

    used_restart_global[i] = 0;
  }

  // nfix_restart_peratom = # of restart entries with peratom info

  int maxsize = 0;

  if (me == 0) utils::sfread(FLERR, &nfix_restart_peratom, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&nfix_restart_peratom, 1, MPI_INT, 0, world);

  // allocate space for each entry

  if (nfix_restart_peratom) {
    id_restart_peratom = new char *[nfix_restart_peratom];
    style_restart_peratom = new char *[nfix_restart_peratom];
    index_restart_peratom = new int[nfix_restart_peratom];
    used_restart_peratom = new int[nfix_restart_peratom];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, maxsize of one atom's data
  // set index = which set of extra data this fix represents

  for (int i = 0; i < nfix_restart_peratom; i++) {
    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    id_restart_peratom[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, id_restart_peratom[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(id_restart_peratom[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    style_restart_peratom[i] = new char[n];
    if (me == 0)
      utils::sfread(FLERR, style_restart_peratom[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(style_restart_peratom[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    maxsize += n;

    index_restart_peratom[i] = i;
    used_restart_peratom[i] = 0;
  }

  return maxsize;
}

/* ----------------------------------------------------------------------
   delete all lists of restart file Fix info
   if flag set, print list of restart file info not assigned to new fixes
------------------------------------------------------------------------- */

void Modify::restart_deallocate(int flag)
{
  if (nfix_restart_global) {
    if (flag && comm->me == 0) {
      int i;
      for (i = 0; i < nfix_restart_global; i++)
        if (used_restart_global[i] == 0) break;
      if (i == nfix_restart_global) {
        utils::logmesg(lmp, "All restart file global fix info was re-assigned\n");
      } else {
        utils::logmesg(lmp, "Unused restart file global fix info:\n");
        for (i = 0; i < nfix_restart_global; i++) {
          if (used_restart_global[i]) continue;
          utils::logmesg(lmp, "  fix style: {}, fix ID: {}\n", style_restart_global[i],
                         id_restart_global[i]);
        }
      }
    }

    for (int i = 0; i < nfix_restart_global; i++) {
      delete[] id_restart_global[i];
      delete[] style_restart_global[i];
      delete[] state_restart_global[i];
    }
    delete[] id_restart_global;
    delete[] style_restart_global;
    delete[] state_restart_global;
    delete[] used_restart_global;
  }

  if (nfix_restart_peratom) {
    if (flag && comm->me == 0) {
      int i;
      for (i = 0; i < nfix_restart_peratom; i++)
        if (used_restart_peratom[i] == 0) break;
      if (i == nfix_restart_peratom) {
        utils::logmesg(lmp, "All restart file peratom fix info was re-assigned\n");
      } else {
        utils::logmesg(lmp, "Unused restart file peratom fix info:\n");
        for (i = 0; i < nfix_restart_peratom; i++) {
          if (used_restart_peratom[i]) continue;
          utils::logmesg(lmp, "  fix style: {}, fix ID: {}\n", style_restart_peratom[i],
                         id_restart_peratom[i]);
        }
      }
    }

    for (int i = 0; i < nfix_restart_peratom; i++) {
      delete[] id_restart_peratom[i];
      delete[] style_restart_peratom[i];
    }
    delete[] id_restart_peratom;
    delete[] style_restart_peratom;
    delete[] index_restart_peratom;
    delete[] used_restart_peratom;
  }

  nfix_restart_global = nfix_restart_peratom = 0;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete[] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for end_of_step fixes
   also create end_of_step_every[]
------------------------------------------------------------------------- */

void Modify::list_init_end_of_step(int mask, int &n, int *&list)
{
  delete[] list;
  delete[] end_of_step_every;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    }
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes that compute reservoir coupling energy
   only added to list if fix has ecouple_flag set
------------------------------------------------------------------------- */

void Modify::list_init_energy_couple(int &n, int *&list)
{
  delete[] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->ecouple_flag) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->ecouple_flag) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes that compute global energy
   only added to list if fix has energy_global_flag and thermo_energy set
------------------------------------------------------------------------- */

void Modify::list_init_energy_global(int &n, int *&list)
{
  delete[] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_global_flag && fix[i]->thermo_energy) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_global_flag && fix[i]->thermo_energy) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes that compute peratom energy
   only added to list if fix has energy_peratom_flag and thermo_energy set
------------------------------------------------------------------------- */

void Modify::list_init_energy_atom(int &n, int *&list)
{
  delete[] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_peratom_flag && fix[i]->thermo_energy) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->energy_peratom_flag && fix[i]->thermo_energy) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for fix GROUP
   are invoked first in post_force() or post_force_respa()
------------------------------------------------------------------------- */

void Modify::list_init_post_force_group(int &n, int *&list)
{
  delete[] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (strcmp(fix[i]->style, "GROUP") == 0) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of compute indices for computes which store invocation times
------------------------------------------------------------------------- */

void Modify::list_init_compute()
{
  delete[] list_timeflag;

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) n_timeflag++;
  list_timeflag = new int[n_timeflag];

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) list_timeflag[n_timeflag++] = i;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory from all fixes
------------------------------------------------------------------------- */

double Modify::memory_usage()
{
  double bytes = 0;
  for (int i = 0; i < nfix; i++) bytes += fix[i]->memory_usage();
  for (int i = 0; i < ncompute; i++) bytes += compute[i]->memory_usage();
  return bytes;
}
