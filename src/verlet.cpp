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
#include "verlet.h"
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
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void Verlet::init()
{
  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning("No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_compute via sum over ghost atoms

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // set flags for what arrays to clear in force_clear()
  // need to clear torques if array exists

  torqueflag = 0;
  if (atom->torque_flag) torqueflag = 1;

  // orthogonal vs triclinic simulation box

  triclinic = domain->triclinic;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void Verlet::setup()
{
  if (comm->me == 0 && screen) fprintf(screen,"Setting up run ...\n");

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
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag,vflag);
  }

  if (force->newton) comm->reverse_communicate();

  modify->setup(vflag);
  output->setup(1);
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void Verlet::setup_minimal(int flag)
{
  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    neighbor->build();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag,vflag);
  }

  if (force->newton) comm->reverse_communicate();

  modify->setup(vflag);
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void Verlet::run(int n)
{
  int nflag,ntimestep;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    modify->initial_integrate(vflag);
    if (n_post_integrate) modify->post_integrate();

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->communicate();
      timer->stamp(TIME_COMM);
    } else {
      if (n_pre_exchange) modify->pre_exchange();
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      if (domain->box_change) {
	domain->reset_box();
	comm->setup();
	if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();
      comm->exchange();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      timer->stamp(TIME_COMM);
      if (n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
      timer->stamp(TIME_NEIGHBOR);
    }

    // force computations

    force_clear();
    if (n_pre_force) modify->pre_force(vflag);

    timer->stamp();

    if (force->pair) {
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }

    if (atom->molecular) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }

    if (force->kspace) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    // reverse communication of forces

    if (force->newton) {
      comm->reverse_communicate();
      timer->stamp(TIME_COMM);
    }

    // force modifications, final time integration, diagnostics

    if (n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (n_end_of_step) modify->end_of_step();

    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void Verlet::force_clear()
{
  int i;

  // clear force on all particles
  // if either newton flag is set, also include ghosts

  if (neighbor->includegroup == 0) {
    int nall;
    if (force->newton) nall = atom->nlocal + atom->nghost;
    else nall = atom->nlocal;

    double **f = atom->f;
    for (i = 0; i < nall; i++) {
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
    }
    
    if (torqueflag) {
      double **torque = atom->torque;
      for (i = 0; i < nall; i++) {
	torque[i][0] = 0.0;
	torque[i][1] = 0.0;
	torque[i][2] = 0.0;
      }
    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    int nall = atom->nfirst;

    double **f = atom->f;
    for (i = 0; i < nall; i++) {
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
    }
    
    if (torqueflag) {
      double **torque = atom->torque;
      for (i = 0; i < nall; i++) {
	torque[i][0] = 0.0;
	torque[i][1] = 0.0;
	torque[i][2] = 0.0;
      }
    }

    if (force->newton) {
      nall = atom->nlocal + atom->nghost;

      for (i = atom->nlocal; i < nall; i++) {
	f[i][0] = 0.0;
	f[i][1] = 0.0;
	f[i][2] = 0.0;
      }
    
      if (torqueflag) {
	double **torque = atom->torque;
	for (i = atom->nlocal; i < nall; i++) {
	  torque[i][0] = 0.0;
	  torque[i][1] = 0.0;
	  torque[i][2] = 0.0;
	}
      }
    }
  }
}
