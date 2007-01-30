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

  if (modify->nfix == 0)
    error->warning("No fixes defined, atoms won't move");

  // set flags for how virial should be computed when needed
  // pressure_flag is 1 if NPT,NPH
  // virial_every is how virial should be computed every timestep
  //   0 = not computed, 1 = computed explicity by pair, 
  //   2 = computed implicitly by pair (via summation over ghost atoms)
  // virial_thermo is how virial should be computed on thermo timesteps
  //   1 = computed explicity by pair, 2 = computed implicitly by pair

  int pressure_flag = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"npt") == 0) pressure_flag = 1;
    if (strcmp(modify->fix[i]->style,"nph") == 0) pressure_flag = 1;
  }

  if (pressure_flag && force->newton_pair) virial_every = 2;
  else if (pressure_flag) virial_every = 1;
  else virial_every = 0;

  if (force->newton_pair) virial_thermo = 2;
  else virial_thermo = 1;

  // set flags for what arrays to clear in force_clear()
  // need to clear torques if atom_style is dipole
  // need to clear phia if atom_style is granular
  // don't need to clear f_pair if atom_style is only granular (no virial)

  torqueflag = 0;
  if (atom->check_style("dipole")) torqueflag = 1;
  granflag = 0;
  if (atom->check_style("granular")) granflag = 1;
  pairflag = 1;
  if (strcmp(atom->atom_style,"granular") == 0) pairflag = 0;

  // local versions of Update quantities

  maxpair = update->maxpair;
  f_pair = update->f_pair;
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

  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  neighbor->build();
  neighbor->ncalls = 0;

  // compute all forces

  int eflag = 1;
  int vflag = virial_thermo;
  force_clear(vflag);
  
  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->pair) force->pair->compute(eflag,vflag);

  if (force->kspace) {
    force->kspace->setup();
    force->kspace->compute(eflag,vflag);
  }

  if (force->newton) comm->reverse_communicate();

  modify->setup();
  output->setup(1);
}

/* ----------------------------------------------------------------------
   iterate for n steps
------------------------------------------------------------------------- */

void Verlet::iterate(int n)
{
  int eflag,vflag,nflag;

  for (int i = 0; i < n; i++) {

    update->ntimestep++;

    modify->initial_integrate();

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->communicate();
      timer->stamp(TIME_COMM);
    } else {
      if (modify->n_pre_exchange) modify->pre_exchange();
      domain->pbc();
      if (domain->box_change) {
	domain->reset_box();
	comm->setup();
	if (neighbor->style) neighbor->setup_bins();
      }
      timer->stamp();
      comm->exchange();
      comm->borders();
      timer->stamp(TIME_COMM);
      if (modify->n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
      timer->stamp(TIME_NEIGHBOR);
    }

    eflag = 0;
    vflag = virial_every;
    if (output->next_thermo == update->ntimestep) {
      eflag = 1;
      vflag = virial_thermo;
    }
    force_clear(vflag);

    timer->stamp();
    if (atom->molecular) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
      timer->stamp(TIME_BOND);
    }

    if (force->pair) {
      force->pair->compute(eflag,vflag);
      timer->stamp(TIME_PAIR);
    }

    if (force->kspace) {
      force->kspace->compute(eflag,vflag);
      timer->stamp(TIME_KSPACE);
    }

    if (force->newton) {
      comm->reverse_communicate();
      timer->stamp(TIME_COMM);
    }

    if (modify->n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (modify->n_end_of_step) modify->end_of_step();

    if (output->next == update->ntimestep) {
      timer->stamp();
      output->write(update->ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   setup and clear other arrays as needed
------------------------------------------------------------------------- */

void Verlet::force_clear(int vflag)
{
  int i;

  // clear global force array
  // nall includes ghosts only if either newton flag is set

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

  if (granflag) {
    double **phia = atom->phia;
    for (i = 0; i < nall; i++) {
      phia[i][0] = 0.0;
      phia[i][1] = 0.0;
      phia[i][2] = 0.0;
    }
  }

  // clear f_pair array if using it this timestep to compute virial

  if (vflag == 2 && pairflag) {
    if (atom->nmax > maxpair) {
      maxpair = atom->nmax;
      memory->destroy_2d_double_array(f_pair);
      f_pair = memory->create_2d_double_array(maxpair,3,"verlet:f_pair");
      update->maxpair = maxpair;
      update->f_pair = f_pair;
    }
    for (i = 0; i < nall; i++) {
      f_pair[i][0] = 0.0;
      f_pair[i][1] = 0.0;
      f_pair[i][2] = 0.0;
    }
  }
}
