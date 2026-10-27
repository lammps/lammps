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

#include "accelerator_kokkos.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"    // IWYU pragma: keep
#include "domain.h"
#include "error.h"
#include "fix.h"    // IWYU pragma: keep
#include "group.h"
#include "input.h"
#include "label_map.h"
#include "memory.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double BIG = 1.0e20;

/* ----------------------------------------------------------------------
   process-global registries of fix and compute style factory functions.
   Shared by all LAMMPS instances and persistent across the "clear" command.
   Built-in styles are registered once by the generated register_*_styles()
   functions; plugins add/override entries at runtime.
------------------------------------------------------------------------- */

CreatorRegistry<Modify::FixCreator> &Modify::fix_styles()
{
  static CreatorRegistry<Modify::FixCreator> registry;
  return registry;
}

CreatorRegistry<Modify::ComputeCreator> &Modify::compute_styles()
{
  static CreatorRegistry<Modify::ComputeCreator> registry;
  return registry;
}

/* ---------------------------------------------------------------------- */

Modify::Modify(LAMMPS *lmp) : Pointers(lmp)
{
  n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = n_post_neighbor = 0;
  n_pre_force = n_pre_reverse = n_post_force_any = 0;
  n_final_integrate = n_end_of_step = 0;
  n_energy_couple = n_energy_global = n_energy_atom = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa_any = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_neighbor = n_min_post_neighbor = 0;
  n_min_pre_force = n_min_pre_reverse = n_min_post_force = n_min_energy = 0;
  n_post_force = n_post_force_group = n_post_force_respa = 0;

  nfix = 0;
  fix = nullptr;
  fmask = nullptr;
  ncompute = 0;
  compute = nullptr;

  restart_pbc_any = 0;
  nfix_restart_global = 0;
  id_restart_global = style_restart_global = nullptr;
  state_restart_global = nullptr;
  used_restart_global = nullptr;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = nullptr;
  index_restart_peratom = used_restart_peratom = nullptr;
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes
  // do it via delete_fix() so callbacks in Atom are also updated correctly

  while (!fixes.empty()) {
    if (fixes[0])
      delete_fix(fixes[0]);
    else
      erase_fix_slot(0);
  }

  // delete all computes
  // do it via delete_compute() for clean deletion of computes that have created other computes

  while (!computes.empty()) delete_compute(computes[0]);

  restart_deallocate(0);
}

/* ----------------------------------------------------------------------
   refresh the temporary raw views of the lists after they were changed
------------------------------------------------------------------------- */

void Modify::refresh_views()
{
  nfix = (int) fixes.size();
  fix = fixes.data();
  fmask = fixmask.data();
  ncompute = (int) computes.size();
  compute = computes.data();
}

/* ----------------------------------------------------------------------
   count fixes whose mask has any of the bits in mask set
------------------------------------------------------------------------- */

int Modify::count_fix_mask(int mask) const
{
  int n = 0;
  for (const auto &m : callmask)
    if (m & mask) n++;
  return n;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate(1);

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes
  // perform check whether extscalar, extvector, and extarray have been
  //   set when scalar_flag, vector_flag, or array_flag are true.
  // loop re-reads the list size since a compute init() may add computes

  for (std::size_t i = 0; i < computes.size(); i++) {
    computes[i]->init();
    computes[i]->init_flags();
  }
  addstep_compute_all(update->ntimestep);

  // init each fix
  // should not need to come before compute init
  //   used to b/c temperature computes called fix->dof() in their init,
  //   and fix rigid required its own init before its dof() could be called,
  //   but computes now do their DOF in setup()
  // perform check whether extscalar, extvector, and extarray have been
  //   set when scalar_flag, vector_flag, or array_flag are true.
  // loop re-reads the list size since a fix init() may add fixes

  for (std::size_t i = 0; i < fixes.size(); i++) {
    fixes[i]->init();
    fixes[i]->init_flags();
  }

  // set global flag if any fix has its restart_pbc flag set

  restart_pbc_any = 0;
  for (const auto &ifix : fixes)
    if (ifix->restart_pbc) restart_pbc_any = 1;

  // take the snapshot of the fix masks used to invoke fixes until the next
  //   init() and count fixes to call at each stage of run
  // needs to happen after init() of computes and fixes
  //   b/c a compute::init() can delete a fix, e.g. compute chunk/atom
  //   and a fix::init() can add fixes, e.g. fix reaxff/species

  callmask = fixmask;

  n_initial_integrate = count_fix_mask(INITIAL_INTEGRATE);
  n_post_integrate = count_fix_mask(POST_INTEGRATE);
  n_pre_exchange = count_fix_mask(PRE_EXCHANGE);
  n_pre_neighbor = count_fix_mask(PRE_NEIGHBOR);
  n_post_neighbor = count_fix_mask(POST_NEIGHBOR);
  n_pre_force = count_fix_mask(PRE_FORCE);
  n_pre_reverse = count_fix_mask(PRE_REVERSE);
  n_post_force = count_fix_mask(POST_FORCE);
  n_post_force_group = count_fix_mask(GROUP_FIX);
  n_final_integrate = count_fix_mask(FINAL_INTEGRATE);
  n_end_of_step = count_fix_mask(END_OF_STEP);

  n_initial_integrate_respa = count_fix_mask(INITIAL_INTEGRATE_RESPA);
  n_post_integrate_respa = count_fix_mask(POST_INTEGRATE_RESPA);
  n_post_force_respa = count_fix_mask(POST_FORCE_RESPA);
  n_pre_force_respa = count_fix_mask(PRE_FORCE_RESPA);
  n_final_integrate_respa = count_fix_mask(FINAL_INTEGRATE_RESPA);

  n_min_pre_exchange = count_fix_mask(MIN_PRE_EXCHANGE);
  n_min_pre_neighbor = count_fix_mask(MIN_PRE_NEIGHBOR);
  n_min_post_neighbor = count_fix_mask(MIN_POST_NEIGHBOR);
  n_min_pre_force = count_fix_mask(MIN_PRE_FORCE);
  n_min_pre_reverse = count_fix_mask(MIN_PRE_REVERSE);
  n_min_post_force = count_fix_mask(MIN_POST_FORCE);
  n_min_energy = count_fix_mask(MIN_ENERGY);

  // two post_force_any counters used by integrators add in post_force_group

  n_post_force_any = n_post_force + n_post_force_group;
  n_post_force_respa_any = n_post_force_respa + n_post_force_group;

  // count fixes that contribute to the energy sums (selected by flags, not mask)

  n_energy_couple = n_energy_global = n_energy_atom = 0;
  for (const auto &ifix : fixes) {
    if (ifix->ecouple_flag) n_energy_couple++;
    if (ifix->energy_global_flag && ifix->thermo_energy) n_energy_global++;
    if (ifix->energy_peratom_flag && ifix->thermo_energy) n_energy_atom++;
  }

  // error if any fix or compute is using a dynamic group when not allowed

  for (const auto &ifix : fixes)
    if (!ifix->dynamic_group_allow && group->dynamic[ifix->igroup])
      error->all(FLERR, Error::NOLASTLINE, "Fix {} does not allow use with a dynamic group",
                 ifix->style);

  for (const auto &icompute : computes)
    if (!icompute->dynamic_group_allow && group->dynamic[icompute->igroup])
      error->all(FLERR, Error::NOLASTLINE, "Compute {} does not allow use with a dynamic group",
                 icompute->style);

  // warn if any particle is time integrated more than once

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int *flag = new int[nlocal];
  for (int i = 0; i < nlocal; i++) flag[i] = 0;

  for (const auto &ifix : fixes) {
    if (ifix->time_integrate == 0) continue;
    const int groupbit = ifix->groupbit;
    for (int j = 0; j < nlocal; j++)
      if (mask[j] & groupbit) flag[j]++;
  }

  int check = 0;
  for (int i = 0; i < nlocal; i++)
    if (flag[i] > 1) check = 1;

  delete[] flag;

  int checkall;
  MPI_Allreduce(&check, &checkall, 1, MPI_INT, MPI_SUM, world);
  if (comm->me == 0 && checkall)
    error->warning(FLERR,
                   "One or more atoms are time integrated more than once" + utils::errorurl(32));
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

  for_each_fix(GROUP_FIX, [vflag](Fix *ifix) {
    ifix->setup(vflag);
  });

  for (std::size_t i = 0; i < computes.size(); i++) computes[i]->setup();

  if (update->whichflag == 1)
    for (std::size_t i = 0; i < fixes.size(); i++) fixes[i]->setup(vflag);
  else if (update->whichflag == 2)
    for (std::size_t i = 0; i < fixes.size(); i++) fixes[i]->min_setup(vflag);

  // runtime check for type label self-consistency

  if (atom->labelmapflag && atom->lmap->checkflag) atom->lmap->check_labels();
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void Modify::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    for_each_fix(PRE_EXCHANGE, [](Fix *ifix) {
      ifix->setup_pre_exchange();
    });
  else if (update->whichflag == 2)
    for_each_fix(MIN_PRE_EXCHANGE, [](Fix *ifix) {
      ifix->setup_pre_exchange();
    });
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    for_each_fix(PRE_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_pre_neighbor();
    });
  else if (update->whichflag == 2)
    for_each_fix(MIN_PRE_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_pre_neighbor();
    });
}

/* ----------------------------------------------------------------------
   setup post_neighbor call, only for fixes that define post_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_post_neighbor()
{
  if (update->whichflag == 1)
    for_each_fix(POST_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_post_neighbor();
    });
  else if (update->whichflag == 2)
    for_each_fix(MIN_POST_NEIGHBOR, [](Fix *ifix) {
      ifix->setup_post_neighbor();
    });
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_force(int vflag)
{
  if (update->whichflag == 1)
    for_each_fix(PRE_FORCE, [vflag](Fix *ifix) {
      ifix->setup_pre_force(vflag);
    });
  else if (update->whichflag == 2)
    for_each_fix(MIN_PRE_FORCE, [vflag](Fix *ifix) {
      ifix->setup_pre_force(vflag);
    });
}

/* ----------------------------------------------------------------------
   setup pre_reverse call, only for fixes that define pre_reverse
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_reverse(int eflag, int vflag)
{
  if (update->whichflag == 1)
    for_each_fix(PRE_REVERSE, [=](Fix *ifix) {
      ifix->setup_pre_reverse(eflag, vflag);
    });
  else if (update->whichflag == 2)
    for_each_fix(MIN_PRE_REVERSE, [=](Fix *ifix) {
      ifix->setup_pre_reverse(eflag, vflag);
    });
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate(int vflag)
{
  for_each_fix(INITIAL_INTEGRATE, [vflag](Fix *ifix) {
    ifix->initial_integrate(vflag);
  });
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate()
{
  for_each_fix(POST_INTEGRATE, [](Fix *ifix) {
    ifix->post_integrate();
  });
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_exchange()
{
  for_each_fix(PRE_EXCHANGE, [](Fix *ifix) {
    ifix->pre_exchange();
  });
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_neighbor()
{
  for_each_fix(PRE_NEIGHBOR, [](Fix *ifix) {
    ifix->pre_neighbor();
  });
}

/* ----------------------------------------------------------------------
   post_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_neighbor()
{
  for_each_fix(POST_NEIGHBOR, [](Fix *ifix) {
    ifix->post_neighbor();
  });
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force(int vflag)
{
  for_each_fix(PRE_FORCE, [vflag](Fix *ifix) {
    ifix->pre_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   pre_reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_reverse(int eflag, int vflag)
{
  for_each_fix(PRE_REVERSE, [=](Fix *ifix) {
    ifix->pre_reverse(eflag, vflag);
  });
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
     they are not in n_post_force count
------------------------------------------------------------------------- */

void Modify::post_force(int vflag)
{
  if (n_post_force_group)
    for_each_fix(GROUP_FIX, [vflag](Fix *ifix) {
      ifix->post_force(vflag);
    });
  if (n_post_force)
    for_each_fix(POST_FORCE, [vflag](Fix *ifix) {
      ifix->post_force(vflag);
    });
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  for_each_fix(FINAL_INTEGRATE, [](Fix *ifix) {
    ifix->final_integrate();
  });
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
  const bigint ntimestep = update->ntimestep;
  for_each_fix(END_OF_STEP, [ntimestep](Fix *ifix) {
    if (ntimestep % ifix->nevery == 0) ifix->end_of_step();
  });
}

/* ----------------------------------------------------------------------
   coupling energy call, only for relevant fixes
   each thermostat fix returns this via compute_scalar()
   ecouple = cumulative energy added to reservoir by thermostatting
------------------------------------------------------------------------- */

double Modify::energy_couple()
{
  double energy = 0.0;
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->ecouple_flag) energy += fixes[i]->compute_scalar();
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
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->energy_global_flag && fixes[i]->thermo_energy)
      energy += fixes[i]->compute_scalar();
  return energy;
}

/* ----------------------------------------------------------------------
   peratom energy call, only for relevant fixes
   called by compute pe/atom
------------------------------------------------------------------------- */

void Modify::energy_atom(int nlocal, double *energy)
{
  for (const auto &ifix : fixes) {
    if (!ifix->energy_peratom_flag || !ifix->thermo_energy) continue;
    const double *eatom = ifix->eatom;
    if (!eatom) continue;
    for (int j = 0; j < nlocal; j++) energy[j] += eatom[j];
  }
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void Modify::post_run()
{
  for (std::size_t i = 0; i < fixes.size(); i++) fixes[i]->post_run();
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
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->create_attribute) fixes[i]->set_arrays(n);
  for (std::size_t i = 0; i < computes.size(); i++)
    if (computes[i]->create_attribute) computes[i]->set_arrays(n);
  input->variable->set_arrays(n);
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::setup_pre_force_respa(int vflag, int ilevel)
{
  for_each_fix(PRE_FORCE_RESPA, [=](Fix *ifix) {
    ifix->setup_pre_force_respa(vflag, ilevel);
  });
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  for_each_fix(INITIAL_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->initial_integrate_respa(vflag, ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate_respa(int ilevel, int iloop)
{
  for_each_fix(POST_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->post_integrate_respa(ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force_respa(int vflag, int ilevel, int iloop)
{
  for_each_fix(PRE_FORCE_RESPA, [=](Fix *ifix) {
    ifix->pre_force_respa(vflag, ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
   first call any instances of fix GROUP if they exist
------------------------------------------------------------------------- */

void Modify::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (n_post_force_group)
    for_each_fix(GROUP_FIX, [=](Fix *ifix) {
      ifix->post_force_respa(vflag, ilevel, iloop);
    });
  if (n_post_force_respa)
    for_each_fix(POST_FORCE_RESPA, [=](Fix *ifix) {
      ifix->post_force_respa(vflag, ilevel, iloop);
    });
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate_respa(int ilevel, int iloop)
{
  for_each_fix(FINAL_INTEGRATE_RESPA, [=](Fix *ifix) {
    ifix->final_integrate_respa(ilevel, iloop);
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_exchange()
{
  for_each_fix(MIN_PRE_EXCHANGE, [](Fix *ifix) {
    ifix->min_pre_exchange();
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_neighbor()
{
  for_each_fix(MIN_PRE_NEIGHBOR, [](Fix *ifix) {
    ifix->min_pre_neighbor();
  });
}

/* ----------------------------------------------------------------------
   minimizer post-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_neighbor()
{
  for_each_fix(MIN_POST_NEIGHBOR, [](Fix *ifix) {
    ifix->min_post_neighbor();
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_force(int vflag)
{
  for_each_fix(MIN_PRE_FORCE, [vflag](Fix *ifix) {
    ifix->min_pre_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer pre-reverse call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_reverse(int eflag, int vflag)
{
  for_each_fix(MIN_PRE_REVERSE, [=](Fix *ifix) {
    ifix->min_pre_reverse(eflag, vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_force(int vflag)
{
  for_each_fix(MIN_POST_FORCE, [vflag](Fix *ifix) {
    ifix->min_post_force(vflag);
  });
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double Modify::min_energy(double *fextra)
{
  int index = 0;
  double eng = 0.0;
  for_each_fix(MIN_ENERGY, [&](Fix *ifix) {
    eng += ifix->min_energy(&fextra[index]);
    index += ifix->min_dof();
  });
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_store()
{
  for_each_fix(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_store();
  });
}

/* ----------------------------------------------------------------------
   manage state of extra minimizer dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_clearstore()
{
  for_each_fix(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_clearstore();
  });
}

void Modify::min_pushstore()
{
  for_each_fix(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_pushstore();
  });
}

void Modify::min_popstore()
{
  for_each_fix(MIN_ENERGY, [](Fix *ifix) {
    ifix->min_popstore();
  });
}

/* ----------------------------------------------------------------------
   displace extra minimizer dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_step(double alpha, double *hextra)
{
  int index = 0;
  for_each_fix(MIN_ENERGY, [&](Fix *ifix) {
    ifix->min_step(alpha, &hextra[index]);
    index += ifix->min_dof();
  });
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double Modify::max_alpha(double *hextra)
{
  int index = 0;
  double alpha = BIG;
  for_each_fix(MIN_ENERGY, [&](Fix *ifix) {
    double alpha_one = ifix->max_alpha(&hextra[index]);
    alpha = MIN(alpha, alpha_one);
    index += ifix->min_dof();
  });
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra minimizer dof, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_dof()
{
  int ndof = 0;
  for_each_fix(MIN_ENERGY, [&](Fix *ifix) {
    ndof += ifix->min_dof();
  });
  return ndof;
}

/* ----------------------------------------------------------------------
   reset minimizer reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_reset_ref()
{
  int itmpall = 0;
  for_each_fix(MIN_ENERGY, [&](Fix *ifix) {
    if (ifix->min_reset_ref()) itmpall = 1;
  });
  return itmpall;
}

/* ----------------------------------------------------------------------
   reset grids for any Fix or Compute that uses distributed grids
   called by load balancer when proc sub-domains change
------------------------------------------------------------------------- */

void Modify::reset_grid()
{
  for (std::size_t i = 0; i < fixes.size(); i++)
    if (fixes[i]->pergrid_flag) fixes[i]->reset_grid();
  for (std::size_t i = 0; i < computes.size(); i++)
    if (computes[i]->pergrid_flag) computes[i]->reset_grid();
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
    {"GPU", "OMP", "INTEL", "property/atom", "cmap", "cmap3", "rx", "deprecated", "STORE/KIM",
     "amoeba/pitorsion", "amoeba/bitorsion", "DUMMY", nullptr};
  // clang-format on

  if (domain->box_exist == 0) {
    int m;
    for (m = 0; exceptions[m] != nullptr; m++)
      if (strcmp(arg[2], exceptions[m]) == 0) break;
    if (exceptions[m] == nullptr)
      error->all(FLERR, 2, "Fix {} command before simulation box is defined", arg[2]);
  }

  // check group ID

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR, 1, "Could not find fix group ID {}", arg[1]);

  // if fix ID exists:
  //   create the new fix in the same slot of the list, so its position
  //     relative to the other fixes (= order of invocation) is kept
  //   error if new style does not match old style
  //   warn if new group != old group
  //   delete old fix, but do not call atom->update_callback(),
  //     since the new fix takes the same slot and no other fix moves
  //   leave a null pointer in the slot while the new fix is constructed,
  //     e.g. Atom::add_callback() called by the new fix scans the list for it
  // if fix ID does not exist:
  //   append the new fix to the list after it has been constructed

  std::size_t ifix = 0;
  for (; ifix < fixes.size(); ifix++)
    if (fixes[ifix] && (strcmp(arg[0], fixes[ifix]->id) == 0)) break;
  const bool replace = (ifix < fixes.size());

  if (replace) {
    Fix *oldfix = fixes[ifix];
    int match = 0;
    if (strcmp(arg[2], oldfix->style) == 0) match = 1;
    if (!match && trysuffix && lmp->suffix_enable) {
      if (lmp->non_pair_suffix()) {
        std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
        if (estyle == oldfix->style) match = 1;
      }
      if (lmp->suffix2) {
        std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
        if (estyle == oldfix->style) match = 1;
      }
    }
    if (!match) error->all(FLERR, 2, "Replacing a fix, but new style != old style");

    if (oldfix->igroup != igroup && comm->me == 0)
      error->warning(FLERR, "Replacing a fix, but new group != old group");
    delete oldfix;
    fixes[ifix] = nullptr;
    fixmask[ifix] = callmask[ifix] = 0;
    refresh_views();
  }

  // if the constructor of the replacement fix fails with an exception
  // (library interface), remove the empty slot so that the list stays valid

  struct SlotGuard {
    Modify *modify;
    std::size_t index;
    bool armed;
    ~SlotGuard()
    {
      if (armed) modify->erase_fix_slot(index);
    }
  } guard{this, ifix, replace};

  // create the Fix
  // try first with suffix appended

  Fix *newfix = nullptr;

  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
      FixCreator fix_creator = fix_styles().find(estyle);
      if (fix_creator) {
        newfix = fix_creator(lmp, narg, arg);
        delete[] newfix->style;
        newfix->style = utils::strdup(estyle);
      }
    }
    if ((newfix == nullptr) && lmp->suffix2) {
      std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
      FixCreator fix_creator = fix_styles().find(estyle);
      if (fix_creator) {
        newfix = fix_creator(lmp, narg, arg);
        delete[] newfix->style;
        newfix->style = utils::strdup(estyle);
      }
    }
  }

  // fix property/atom keeps its per-atom data in arrays that AtomKokkos owns and
  // grows as Kokkos views.  The plain style would hand those pointers to
  // memory->grow(), i.e. srealloc() on a Kokkos allocation, which corrupts the
  // heap.  Use the KOKKOS version whenever the package is active, even where the
  // suffix machinery above did not apply because suffixes are disabled.

  if ((newfix == nullptr) && lmp->kokkos && (strcmp(arg[2], "property/atom") == 0)) {
    if (FixCreator fix_creator = fix_styles().find("property/atom/kk")) {
      newfix = fix_creator(lmp, narg, arg);
      delete[] newfix->style;
      newfix->style = utils::strdup("property/atom/kk");
    }
  }

  if (newfix == nullptr) {
    if (FixCreator fix_creator = fix_styles().find(arg[2])) newfix = fix_creator(lmp, narg, arg);
  }

  if (newfix == nullptr) error->all(FLERR, 2, utils::check_packages_for_style("fix", arg[2], lmp));

  // store the new fix and its mask in the list
  // the mask is set before post_constructor() so that changes made there
  //   through set_fix_mask() or clear_fix_mask() are kept
  // the new fix is not invoked before the next init() (callmask = 0)

  int mask = newfix->setmask();
  if (strcmp(newfix->style, "GROUP") == 0) mask |= GROUP_FIX;

  if (replace) {
    fixes[ifix] = newfix;
    fixmask[ifix] = mask;
    callmask[ifix] = 0;
  } else {
    fixes.push_back(newfix);
    fixmask.push_back(mask);
    callmask.push_back(0);
  }
  guard.armed = false;
  refresh_views();

  if (newfix->kokkosable && (!lmp->kokkos || !lmp->kokkos->kokkos_exists))
    error->all(FLERR, Error::NOLASTLINE, "Cannot use KOKKOS styles without enabling KOKKOS");

  // post_constructor() can call virtual methods in parent or child
  //   which would otherwise not yet be visible in child class
  // post_constructor() allows new fix to create other fixes
  // it may also delete fixes, so the position of the new fix in the list
  //   may change: only use the pointer to it from here on

  newfix->post_constructor();

  // check if Fix is in restart_global list
  // if yes, pass state info to the Fix so it can reset itself

  for (int i = 0; i < nfix_restart_global; i++)
    if ((strcmp(id_restart_global[i], newfix->id) == 0) &&
        (utils::strip_style_suffix(newfix->style, lmp) == style_restart_global[i])) {
      newfix->restart(state_restart_global[i]);
      used_restart_global[i] = 1;
      newfix->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting global fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       newfix->style, newfix->id);
    }

  // check if Fix is in restart_peratom list
  // if yes, loop over atoms so they can extract info from atom->extra array

  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strcmp(id_restart_peratom[i], newfix->id) == 0 &&
        strcmp(style_restart_peratom[i], newfix->style) == 0) {
      used_restart_peratom[i] = 1;
      for (int j = 0; j < atom->nlocal; j++) newfix->unpack_restart(j, index_restart_peratom[i]);
      newfix->restart_reset = 1;
      if (comm->me == 0)
        utils::logmesg(lmp,
                       "Resetting peratom fix info from restart file:\n"
                       "  fix style: {}, fix ID: {}\n",
                       newfix->style, newfix->id);
    }

  refresh_views();

  // return pointer to fix

  return newfix;
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

Fix *Modify::replace_fix(const std::string &replaceID, int narg, char **arg, int trysuffix)
{
  auto *oldfix = get_fix_by_id(replaceID);
  if (!oldfix)
    error->all(FLERR, Error::NOLASTLINE, "Modify replace_fix ID {} could not be found", replaceID);

  // change ID, igroup, style of fix being replaced to match new fix
  // requires some error checking on arguments for new fix

  if (narg < 3)
    error->all(FLERR, Error::NOLASTLINE, "Not enough arguments for replace_fix invocation");
  if (get_fix_by_id(arg[0])) error->all(FLERR, "Replace_fix ID {} is already in use", arg[0]);

  delete[] oldfix->id;
  oldfix->id = utils::strdup(arg[0]);

  int jgroup = group->find(arg[1]);
  if (jgroup == -1)
    error->all(FLERR, Error::NOLASTLINE, "Could not find replace_fix group ID {}", arg[1]);
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
  return replace_fix(oldfix, args.size(), newarg.data(), trysuffix);
}

/* ----------------------------------------------------------------------
   modify a Fix's parameters
------------------------------------------------------------------------- */

void Modify::modify_fix(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify", error);

  auto *ifix = get_fix_by_id(arg[0]);
  if (!ifix) error->all(FLERR, Error::NOLASTLINE, "Could not find fix_modify ID {}", arg[0]);
  ifix->modify_params(narg - 1, &arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
   Atom class must update indices in its list of callbacks to fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(const std::string &id)
{
  // no more fixes, nothing to do
  if (fixes.empty()) return;

  auto *ifix = get_fix_by_id(id);
  if (!ifix) error->all(FLERR, Error::NOLASTLINE, "Could not find fix ID {} to delete", id);
  delete_fix(ifix);
}

void Modify::delete_fix(Fix *ifix)
{
  if (!ifix) return;

  // locate fix in list; nothing to do if it is not in the list

  std::size_t index = 0;
  while ((index < fixes.size()) && (fixes[index] != ifix)) index++;
  if (index == fixes.size()) return;

  // delete instance while it is still in the list: its destructor looks up
  //   its own position to remove its callbacks in Atom and may delete other
  //   fixes, so the list position is looked up again by pointer value afterwards

  delete ifix;

  if ((index >= fixes.size()) || (fixes[index] != ifix)) {
    index = 0;
    while ((index < fixes.size()) && (fixes[index] != ifix)) index++;
    if (index == fixes.size()) return;
  }
  erase_fix_slot(index);
}

/* ----------------------------------------------------------------------
   remove slot from the fix list and update Atom callbacks to fixes after it
------------------------------------------------------------------------- */

void Modify::erase_fix_slot(std::size_t index)
{
  fixes.erase(fixes.begin() + index);
  fixmask.erase(fixmask.begin() + index);
  callmask.erase(callmask.begin() + index);
  atom->update_callback((int) index);
  refresh_views();
}

/* ----------------------------------------------------------------------
   find a fix by ID
   return index of fix or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_fix(const std::string &id) const
{
  if (id.empty()) return -1;
  for (std::size_t ifix = 0; ifix < fixes.size(); ifix++)
    if (fixes[ifix] && (id == fixes[ifix]->id)) return (int) ifix;
  return -1;
}

/* ----------------------------------------------------------------------
   look up pointer to Fix class by fix-ID
   return null pointer if ID not found
------------------------------------------------------------------------- */

Fix *Modify::get_fix_by_id(const std::string &id) const
{
  if (id.empty()) return nullptr;
  for (const auto &ifix : fixes)
    if (ifix && (id == ifix->id)) return ifix;
  return nullptr;
}

/* ----------------------------------------------------------------------
   look up pointer to fixes by fix style name
   return vector of matching pointers
------------------------------------------------------------------------- */

std::vector<Fix *> Modify::get_fix_by_style(const std::string &style) const
{
  std::vector<Fix *> matches;
  if (style.empty()) return matches;

  for (const auto &ifix : fixes) {
    if (ifix && utils::strmatch(ifix->style, style)) matches.push_back(ifix);
  }
  return matches;
}

/* ----------------------------------------------------------------------
   query and change the mask of stages at which a fix is invoked
   internal bits are not visible and cannot be changed
------------------------------------------------------------------------- */

int Modify::get_fix_mask(Fix *ifix) const
{
  for (std::size_t i = 0; i < fixes.size(); ++i) {
    if (fixes[i] == ifix) return fixmask[i] & PUBLIC_MASK;
  }
  return 0;
}

void Modify::set_fix_mask(Fix *ifix, int flags)
{
  for (std::size_t i = 0; i < fixes.size(); ++i) {
    if (fixes[i] == ifix) fixmask[i] |= (flags & PUBLIC_MASK);
  }
}

void Modify::clear_fix_mask(Fix *ifix, int flags)
{
  for (std::size_t i = 0; i < fixes.size(); ++i) {
    if (fixes[i] == ifix) fixmask[i] &= ~(flags & PUBLIC_MASK);
  }
}

/* ----------------------------------------------------------------------
   check for fix associated with package name in compiled list
   return 1 if found else 0
   used to determine whether LAMMPS was built with
     GPU, INTEL, OPENMP packages, which have their own fixes
------------------------------------------------------------------------- */

int Modify::check_package(const char *package_fix_name)
{
  if (!fix_styles().contains(package_fix_name)) return 0;
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
  for (const auto &ifix : fixes) {
    if (utils::strmatch(ifix->style, "^rigid")) {
      const int *const body = (const int *) ifix->extract("body", dim);
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
  for (const auto &ifix : fixes) {
    if (strncmp("rigid", ifix->style, 5) == 0) {
      const int *const body = (const int *) ifix->extract("body", dim);
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
  for (const auto &ifix : fixes) {
    if (utils::strmatch(ifix->style, "^rigid")) {
      const int *const body = (const int *) ifix->extract("body", dim);
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

  if (get_compute_by_id(arg[0]))
    error->all(FLERR, Error::ARGZERO, "Reuse of compute ID '{}'", arg[0]);

  // create the Compute
  // try first with suffix appended

  Compute *newcompute = nullptr;

  if (trysuffix && lmp->suffix_enable) {
    if (lmp->non_pair_suffix()) {
      std::string estyle = arg[2] + std::string("/") + lmp->non_pair_suffix();
      ComputeCreator compute_creator = compute_styles().find(estyle);
      if (compute_creator) {
        newcompute = compute_creator(lmp, narg, arg);
        delete[] newcompute->style;
        newcompute->style = utils::strdup(estyle);
      }
    }
    if (newcompute == nullptr && lmp->suffix2) {
      std::string estyle = arg[2] + std::string("/") + lmp->suffix2;
      ComputeCreator compute_creator = compute_styles().find(estyle);
      if (compute_creator) {
        newcompute = compute_creator(lmp, narg, arg);
        delete[] newcompute->style;
        newcompute->style = utils::strdup(estyle);
      }
    }
  }

  if (newcompute == nullptr) {
    if (ComputeCreator compute_creator = compute_styles().find(arg[2]))
      newcompute = compute_creator(lmp, narg, arg);
  }

  if (newcompute == nullptr)
    error->all(FLERR, Error::NOLASTLINE, utils::check_packages_for_style("compute", arg[2], lmp));

  // append to list of computes

  computes.push_back(newcompute);
  refresh_views();

  if (newcompute->kokkosable && (!lmp->kokkos || !lmp->kokkos->kokkos_exists))
    error->all(FLERR, Error::NOLASTLINE, "Cannot use KOKKOS styles without enabling KOKKOS");

  // post_constructor() can call virtual methods in parent or child
  //   which would otherwise not yet be visible in child class
  // post_constructor() allows new compute to create other computes
  // the new compute is in the list before, so a recursive call to
  //   add_compute() from post_constructor() will see it

  newcompute->post_constructor();
  refresh_views();

  return newcompute;
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

  auto *icompute = get_compute_by_id(arg[0]);
  if (!icompute)
    error->all(FLERR, Error::NOLASTLINE, "Could not find compute_modify ID {}", arg[0]);
  icompute->modify_params(narg - 1, &arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Compute from list of Computes
------------------------------------------------------------------------- */

void Modify::delete_compute(const std::string &id)
{
  // no more computes, nothing to do
  if (computes.empty()) return;

  auto *icompute = get_compute_by_id(id);
  if (!icompute) error->all(FLERR, Error::NOLASTLINE, "Could not find compute ID {} to delete", id);
  delete_compute(icompute);
}

void Modify::delete_compute(Compute *icompute)
{
  if (!icompute) return;

  // locate compute in list; nothing to do if it is not in the list

  std::size_t index = 0;
  while ((index < computes.size()) && (computes[index] != icompute)) index++;
  if (index == computes.size()) return;

  // delete instance while it is still in the list: its destructor may delete
  //   other computes, so the list position is looked up again by pointer value

  delete icompute;

  if ((index >= computes.size()) || (computes[index] != icompute)) {
    index = 0;
    while ((index < computes.size()) && (computes[index] != icompute)) index++;
    if (index == computes.size()) return;
  }
  computes.erase(computes.begin() + index);
  refresh_views();
}

/* ----------------------------------------------------------------------
   find a compute by ID
   return index of compute or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_compute(const std::string &id) const
{
  if (id.empty()) return -1;
  for (std::size_t icompute = 0; icompute < computes.size(); icompute++)
    if (computes[icompute] && (id == computes[icompute]->id)) return (int) icompute;
  return -1;
}

/* ----------------------------------------------------------------------
   look up pointer to Compute class by compute-ID
   return null pointer if ID not found
------------------------------------------------------------------------- */

Compute *Modify::get_compute_by_id(const std::string &id) const
{
  if (id.empty()) return nullptr;
  for (const auto &icompute : computes)
    if (icompute && (id == icompute->id)) return icompute;
  return nullptr;
}

/* ----------------------------------------------------------------------
   look up pointers to computes by compute style name
   return vector with matching pointers
------------------------------------------------------------------------- */

std::vector<Compute *> Modify::get_compute_by_style(const std::string &style) const
{
  std::vector<Compute *> matches;
  if (style.empty()) return matches;

  for (const auto &icompute : computes) {
    if (icompute && utils::strmatch(icompute->style, style)) matches.push_back(icompute);
  }
  return matches;
}

/* ----------------------------------------------------------------------
   clear invoked flag of all computes
   called everywhere that computes are used, before computes are invoked
   invoked flag used to avoid re-invoking same compute multiple times
   and to flag computes that store invocation times as having been invoked
------------------------------------------------------------------------- */

void Modify::clearstep_compute()
{
  for (const auto &icompute : computes) icompute->invoked_flag = Compute::INVOKED_NONE;
}

/* ----------------------------------------------------------------------
   loop over computes that store invocation times and schedule next invocation
   called everywhere that computes are used, after computes are invoked
   this used to only consider computes listed at init() time whose invoked
     flag was set, but since the flag is never negative that selected every
     listed compute, so now all computes with timeflag set are scheduled
------------------------------------------------------------------------- */

void Modify::addstep_compute(bigint newstep)
{
  addstep_compute_all(newstep);
}

/* ----------------------------------------------------------------------
   loop over all computes
   schedule next invocation for those that store invocation times
   called when not sure what computes will be needed on newstep
------------------------------------------------------------------------- */

void Modify::addstep_compute_all(bigint newstep)
{
  for (std::size_t i = 0; i < computes.size(); i++)
    if (computes[i]->timeflag) computes[i]->addstep(newstep);
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
  for (const auto &ifix : fixes)
    if (ifix->restart_global) count++;

  if (me == 0) fwrite(&count, sizeof(int), 1, fp);

  int n;
  for (const auto &ifix : fixes)
    if (ifix->restart_global) {
      if (me == 0) {
        n = strlen(ifix->id) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(ifix->id, sizeof(char), n, fp);
        auto fix_style = utils::strip_style_suffix(ifix->style, lmp);
        n = fix_style.size() + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(fix_style.c_str(), sizeof(char), n, fp);
      }
      ifix->write_restart(fp);
    }

  count = 0;
  for (const auto &ifix : fixes)
    if (ifix->restart_peratom) count++;

  if (me == 0) fwrite(&count, sizeof(int), 1, fp);

  for (const auto &ifix : fixes)
    if (ifix->restart_peratom) {
      int maxsize_restart = ifix->maxsize_restart();
      if (me == 0) {
        n = strlen(ifix->id) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(ifix->id, sizeof(char), n, fp);
        n = strlen(ifix->style) + 1;
        fwrite(&n, sizeof(int), 1, fp);
        fwrite(ifix->style, sizeof(char), n, fp);
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
  if ((nfix_restart_global < 0) || (nfix_restart_global > 4096))
    error->all(FLERR, "Invalid number of fix entries in restart file");

  // allocate space for each entry

  if (nfix_restart_global > 0) {
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
    if ((n < 1) || (n > 65536)) error->all(FLERR, "Invalid fix info size in restart file");
    id_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, id_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(id_restart_global[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if ((n < 1) || (n > 65536)) error->all(FLERR, "Invalid fix info size in restart file");
    style_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, style_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(style_restart_global[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if ((n < 0) || (n > (1 << 30))) error->all(FLERR, "Invalid fix data size in restart file");
    state_restart_global[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, state_restart_global[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(state_restart_global[i], n, MPI_CHAR, 0, world);

    used_restart_global[i] = 0;
  }

  // nfix_restart_peratom = # of restart entries with peratom info

  bigint maxsize = 0;

  if (me == 0) utils::sfread(FLERR, &nfix_restart_peratom, sizeof(int), 1, fp, nullptr, error);
  MPI_Bcast(&nfix_restart_peratom, 1, MPI_INT, 0, world);
  if ((nfix_restart_peratom < 0) || (nfix_restart_peratom > 4096))
    error->all(FLERR, "Invalid number of fix entries in restart file");

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
    if ((n < 1) || (n > 65536)) error->all(FLERR, "Invalid fix info size in restart file");
    id_restart_peratom[i] = new char[n];
    if (me == 0) utils::sfread(FLERR, id_restart_peratom[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(id_restart_peratom[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if ((n < 1) || (n > 65536)) error->all(FLERR, "Invalid fix info size in restart file");
    style_restart_peratom[i] = new char[n];
    if (me == 0)
      utils::sfread(FLERR, style_restart_peratom[i], sizeof(char), n, fp, nullptr, error);
    MPI_Bcast(style_restart_peratom[i], n, MPI_CHAR, 0, world);

    if (me == 0) utils::sfread(FLERR, &n, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (n < 0) error->all(FLERR, "Invalid per-atom data size in restart file");
    maxsize += n;

    index_restart_peratom[i] = i;
    used_restart_peratom[i] = 0;
  }

  if (maxsize > MAXSMALLINT) error->all(FLERR, "Per-atom fix data in restart file is too large");

  return (int) maxsize;
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
   return # of bytes of allocated memory from all fixes and computes
------------------------------------------------------------------------- */

double Modify::memory_usage()
{
  double bytes = 0;
  for (const auto &ifix : fixes) bytes += ifix->memory_usage();
  for (const auto &icompute : computes) bytes += icompute->memory_usage();
  return bytes;
}
