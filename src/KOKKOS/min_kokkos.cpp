// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "min_kokkos.h"

#include "angle.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "bond.h"
#include "comm.h"
#include "compute.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix_minimize_kokkos.h"
#include "force.h"
#include "improper.h"
#include "kokkos.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "pair.h"
#include "thermo.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MinKokkos::MinKokkos(LAMMPS *lmp) : Min(lmp)
{
  atomKK = (AtomKokkos *) atom;
  fix_minimize_kk = nullptr;
}

/* ---------------------------------------------------------------------- */

MinKokkos::~MinKokkos()
{

}

/* ---------------------------------------------------------------------- */

void MinKokkos::init()
{
  Min::init();

  fix_minimize_kk = (FixMinimizeKokkos*) fix_minimize;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void MinKokkos::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fmt::print(screen,"Setting up {} style minimization ...\n",
               update->minimize_style);
    if (flag) {
      fmt::print(screen,"  Unit style    : {}\n", update->unit_style);
      fmt::print(screen,"  Current step  : {}\n", update->ntimestep);
      timer->print_timeout(screen);
    }
  }
  update->setupflag = 1;

  // setup extra global dof due to fixes
  // cannot be done in init() b/c update init() is before modify init()

  nextra_global = modify->min_dof();
  if (nextra_global) {
    fextra = new double[nextra_global];
    if (comm->me == 0 && screen)
      fprintf(screen,"WARNING: Energy due to %d extra global DOFs will"
              " be included in minimizer energies\n",nextra_global);
  }

  // compute for potential energy

  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Minimization could not find thermo_pe compute");
  pe_compute = modify->compute[id];

  // style-specific setup does two tasks
  // setup extra global dof vectors
  // setup extra per-atom dof vectors due to requests from Pair classes
  // cannot be done in init() b/c update init() is before modify/pair init()

  setup_style();

  // ndoftotal = total dof for entire minimization problem
  // dof for atoms, extra per-atom, extra global

  bigint ndofme = 3 * static_cast<bigint>(atom->nlocal);
  for (int m = 0; m < nextra_atom; m++)
    ndofme += extra_peratom[m]*static_cast<bigint>(atom->nlocal);
  MPI_Allreduce(&ndofme,&ndoftotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
  ndoftotal += nextra_global;

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

  // remove these restriction eventually

  if (searchflag == 0) {
    if (nextra_global)
      error->all(FLERR,
                 "Cannot use a damped dynamics min style with fix box/relax");
    if (nextra_atom)
      error->all(FLERR,
                 "Cannot use a damped dynamics min style with per-atom DOF");
  }

  if (strcmp(update->minimize_style,"hftn") == 0) {
    if (nextra_global)
      error->all(FLERR, "Cannot use hftn min style with fix box/relax");
    if (nextra_atom)
      error->all(FLERR, "Cannot use hftn min style with per-atom DOF");
  }

  // atoms may have migrated in comm->exchange()

  reset_vectors();

  // compute all forces

  force->setup();
  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
    }
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      force->kspace->compute(eflag,vflag);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
    } else force->kspace->compute_dummy(eflag,vflag);
  }

  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  // update per-atom minimization variables stored by pair styles

  if (nextra_atom)
    for (int m = 0; m < nextra_atom; m++)
      requestor[m]->min_xf_get(m);

  lmp->kokkos->auto_sync = 0;
  modify->setup(vflag);
  output->setup(flag);
  lmp->kokkos->auto_sync = 1;
  update->setupflag = 0;

  // stats for initial thermo output

  ecurrent = pe_compute->compute_scalar();
  if (nextra_global) ecurrent += modify->min_energy(fextra);
  if (output->thermo->normflag) ecurrent /= atom->natoms;

  einitial = ecurrent;
  fnorm2_init = sqrt(fnorm_sqr());
  fnorminf_init = sqrt(fnorm_inf());
}

/* ----------------------------------------------------------------------
   setup without output or one-time post-init setup
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void MinKokkos::setup_minimal(int flag)
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

  // atoms may have migrated in comm->exchange()

  reset_vectors();

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
  }
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atomKK->molecular) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
    }
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) {
      atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
      force->kspace->compute(eflag,vflag);
      atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
    } else force->kspace->compute_dummy(eflag,vflag);
  }

  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  // update per-atom minimization variables stored by pair styles

  if (nextra_atom)
    for (int m = 0; m < nextra_atom; m++)
      requestor[m]->min_xf_get(m);

  lmp->kokkos->auto_sync = 0;
  modify->setup(vflag);
  lmp->kokkos->auto_sync = 1;
  update->setupflag = 0;

  // stats for Finish to print

  ecurrent = pe_compute->compute_scalar();
  if (nextra_global) ecurrent += modify->min_energy(fextra);
  if (output->thermo->normflag) ecurrent /= atom->natoms;

  einitial = ecurrent;
  fnorm2_init = sqrt(fnorm_sqr());
  fnorminf_init = sqrt(fnorm_inf());
}

/* ----------------------------------------------------------------------
   perform minimization, calling iterate() for N steps
------------------------------------------------------------------------- */

void MinKokkos::run(int n)
{
  if (nextra_atom)
    error->all(FLERR,"Cannot yet use extra atom DOFs (e.g. AWPMD and EFF packages) "
     "with Kokkos minimize");

  // minimizer iterations

  lmp->kokkos->auto_sync = 0;
  atomKK->sync(Device,ALL_MASK);

  stop_condition = iterate(n);
  stopstr = stopstrings(stop_condition);

  // if early exit from iterate loop:
  // set update->nsteps to niter for Finish stats to print
  // set output->next values to this timestep
  // call energy_force() to insure vflag is set when forces computed
  // output->write does final output for thermo, dump, restart files
  // add ntimestep to all computes that store invocation times
  //   since are hardwiring call to thermo/dumps and computes may not be ready

  if (stop_condition != MAXITER) {
    update->nsteps = niter;

    if (update->restrict_output == 0) {
      for (int idump = 0; idump < output->ndump; idump++)
        output->next_dump[idump] = update->ntimestep;
      output->next_dump_any = update->ntimestep;
      if (output->restart_flag) {
        output->next_restart = update->ntimestep;
        if (output->restart_every_single)
          output->next_restart_single = update->ntimestep;
        if (output->restart_every_double)
          output->next_restart_double = update->ntimestep;
      }
    }
    output->next_thermo = update->ntimestep;

    modify->addstep_compute_all(update->ntimestep);
    ecurrent = energy_force(0);

    atomKK->sync(Host,ALL_MASK);
    output->write(update->ntimestep);
  }

  atomKK->sync(Host,ALL_MASK);
  lmp->kokkos->auto_sync = 1;
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

double MinKokkos::energy_force(int resetflag)
{
  // check for reneighboring
  // always communicate since minimizer moved atoms

  int nflag = neighbor->decide();

  if (nflag == 0) {
    timer->stamp();
    comm->forward_comm();
    timer->stamp(Timer::COMM);
  } else {
    if (modify->n_min_pre_exchange) {
      timer->stamp();
      modify->min_pre_exchange();
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
    if (modify->n_min_pre_neighbor) {
      modify->min_pre_neighbor();
      timer->stamp(Timer::MODIFY);
    }
    neighbor->build(1);
    timer->stamp(Timer::NEIGH);
    if (modify->n_min_post_neighbor) {
      modify->min_post_neighbor();
      timer->stamp(Timer::MODIFY);
    }
  }

  ev_set(update->ntimestep);
  force_clear();

  timer->stamp();

  if (modify->n_min_pre_force) {
    modify->min_pre_force(vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (pair_compute_flag) {
    atomKK->sync(force->pair->execution_space,force->pair->datamask_read);
    force->pair->compute(eflag,vflag);
    atomKK->modified(force->pair->execution_space,force->pair->datamask_modify);
    timer->stamp(Timer::PAIR);
  }

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) {
      atomKK->sync(force->bond->execution_space,force->bond->datamask_read);
      force->bond->compute(eflag,vflag);
      atomKK->modified(force->bond->execution_space,force->bond->datamask_modify);
    }
    if (force->angle) {
      atomKK->sync(force->angle->execution_space,force->angle->datamask_read);
      force->angle->compute(eflag,vflag);
      atomKK->modified(force->angle->execution_space,force->angle->datamask_modify);
    }
    if (force->dihedral) {
      atomKK->sync(force->dihedral->execution_space,force->dihedral->datamask_read);
      force->dihedral->compute(eflag,vflag);
      atomKK->modified(force->dihedral->execution_space,force->dihedral->datamask_modify);
    }
    if (force->improper) {
      atomKK->sync(force->improper->execution_space,force->improper->datamask_read);
      force->improper->compute(eflag,vflag);
      atomKK->modified(force->improper->execution_space,force->improper->datamask_modify);
    }
    timer->stamp(Timer::BOND);
  }

  if (kspace_compute_flag) {
    atomKK->sync(force->kspace->execution_space,force->kspace->datamask_read);
    force->kspace->compute(eflag,vflag);
    atomKK->modified(force->kspace->execution_space,force->kspace->datamask_modify);
    timer->stamp(Timer::KSPACE);
  }

  if (modify->n_min_pre_reverse) {
    modify->min_pre_reverse(eflag,vflag);
    timer->stamp(Timer::MODIFY);
  }

  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }

  // update per-atom minimization variables stored by pair styles

  if (nextra_atom)
    for (int m = 0; m < nextra_atom; m++)
      requestor[m]->min_xf_get(m);

  // fixes that affect minimization

  if (modify->n_min_post_force) {
     timer->stamp();
     modify->min_post_force(vflag);
     timer->stamp(Timer::MODIFY);
  }

  // compute potential energy of system
  // normalize if thermo PE does

  atomKK->sync(pe_compute->execution_space,pe_compute->datamask_read);
  double energy = pe_compute->compute_scalar();
  atomKK->modified(pe_compute->execution_space,pe_compute->datamask_modify);
  if (nextra_global) energy += modify->min_energy(fextra);
  if (output->thermo->normflag) energy /= atom->natoms;

  // if reneighbored, atoms migrated
  // if resetflag = 1, update x0 of atoms crossing PBC
  // reset vectors used by lo-level minimizer

  if (nflag) {
    if (resetflag) fix_minimize_kk->reset_coords();
    reset_vectors();
  }
  return energy;
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void MinKokkos::force_clear()
{
  if (external_force_clear) return;

  // clear global force array
  // if either newton flag is set, also include ghosts

  atomKK->k_f.clear_sync_state(); // ignore host forces/torques since device views
  atomKK->k_torque.clear_sync_state(); // will be cleared below

  int nzero = atom->nlocal;
  if (force->newton) nzero += atom->nghost;

  if (nzero) {
    // local variables for lambda capture

    auto l_f = atomKK->k_f.d_view;
    auto l_torque = atomKK->k_torque.d_view;
    auto l_torqueflag = torqueflag;

    Kokkos::parallel_for(nzero, LAMMPS_LAMBDA(int i) {
      l_f(i,0) = 0.0;
      l_f(i,1) = 0.0;
      l_f(i,2) = 0.0;
      if (l_torqueflag) {
        l_torque(i,0) = 0.0;
        l_torque(i,1) = 0.0;
        l_torque(i,2) = 0.0;
      }
    });
  }
  atomKK->modified(Device,F_MASK);
}

/* ----------------------------------------------------------------------
   compute and return ||force||_2^2
------------------------------------------------------------------------- */

double MinKokkos::fnorm_sqr()
{

  double local_norm2_sqr = 0.0;
  {
    // local variables for lambda capture

    auto l_fvec = fvec;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(int i, double& local_norm2_sqr) {
      local_norm2_sqr += l_fvec[i]*l_fvec[i];
    },local_norm2_sqr);
  }

  double norm2_sqr = 0.0;
  MPI_Allreduce(&local_norm2_sqr,&norm2_sqr,1,MPI_DOUBLE,MPI_SUM,world);

  if (nextra_global)
    for (int i = 0; i < nextra_global; i++)
      norm2_sqr += fextra[i]*fextra[i];

  return norm2_sqr;
}

/* ----------------------------------------------------------------------
   compute and return ||force||_inf
------------------------------------------------------------------------- */

double MinKokkos::fnorm_inf()
{

  double local_norm_inf = 0.0;
  {
    // local variables for lambda capture

    auto l_fvec = fvec;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(int i, double& local_norm_inf) {
      local_norm_inf = MAX(l_fvec[i]*l_fvec[i],local_norm_inf);
    },Kokkos::Max<double>(local_norm_inf));
  }

  double norm_inf = 0.0;
  MPI_Allreduce(&local_norm_inf,&norm_inf,1,MPI_DOUBLE,MPI_MAX,world);

  if (nextra_global)
    for (int i = 0; i < nextra_global; i++)
      norm_inf = MAX(fextra[i]*fextra[i],norm_inf);

  return norm_inf;
}

/* ----------------------------------------------------------------------
   compute and return ||force||_max (inf norm per-vector)
------------------------------------------------------------------------- */

double MinKokkos::fnorm_max()
{

  double local_norm_max = 0.0;
  {
    // local variables for lambda capture

    auto l_fvec = fvec;

    Kokkos::parallel_reduce(nvec, LAMMPS_LAMBDA(int i, double& local_norm_max) {
      double fdotf = l_fvec[i]*l_fvec[i]+l_fvec[i+1]*l_fvec[i+1]+l_fvec[i+2]*l_fvec[i+2];
      local_norm_max = MAX(fdotf,local_norm_max);
    },Kokkos::Max<double>(local_norm_max));
  }

  double norm_max = 0.0;
  MPI_Allreduce(&local_norm_max,&norm_max,1,MPI_DOUBLE,MPI_MAX,world);

  if (nextra_global) {
    for (int i = 0; i < nextra_global; i+=3) {
      double fdotf = fextra[i]*fextra[i];
      norm_max = MAX(fdotf,norm_max);
    }
  }
  return norm_max;
}
