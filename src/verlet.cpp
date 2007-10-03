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

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Verlet::Verlet(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg)
{
  fix_virial_every = NULL;
  next_fix_virial = NULL;
}

/* ---------------------------------------------------------------------- */

Verlet::~Verlet()
{
  delete [] fix_virial_every;
  delete [] next_fix_virial;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void Verlet::init()
{
  // warn if no fixes

  if (modify->nfix == 0)
    error->warning("No fixes defined, atoms won't move");

  // setup virial computations for timestepping
  // virial_style = 1 if computed explicitly by pair
  //                2 if computed implicitly by pair (sum over ghost atoms)
  // virial_every = 1 if computed every timestep (NPT,NPH)
  // fix arrays store info on fixes that need virial computed occasionally

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  virial_every = 0;
  nfix_virial = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->pressure_every == 1) virial_every = 1;
    else if (modify->fix[i]->pressure_every > 1) nfix_virial++;

  if (nfix_virial) {
    delete [] fix_virial_every;
    delete [] next_fix_virial;
    fix_virial_every = new int[nfix_virial];
    next_fix_virial = new int[nfix_virial];
    nfix_virial = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->pressure_every > 1)
	fix_virial_every[nfix_virial++] = modify->fix[i]->pressure_every;
  }

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

  int eflag = 1;
  int vflag = virial_style;
  force_clear(vflag);
  
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

  modify->setup();
  output->setup(1);

  // setup virial computations for timestepping

  int ntimestep = update->ntimestep;
  next_virial = 0;
  if (virial_every) next_virial = ntimestep + 1;
  else {
    for (int ivirial = 0; ivirial < nfix_virial; ivirial++) {
      next_fix_virial[ivirial] = 
	(ntimestep/fix_virial_every[ivirial])*fix_virial_every[ivirial] + 
	fix_virial_every[ivirial];
      if (ivirial) next_virial = MIN(next_virial,next_fix_virial[ivirial]);
      else next_virial = next_fix_virial[0];
    }
  }
}

/* ----------------------------------------------------------------------
   iterate for n steps
------------------------------------------------------------------------- */

void Verlet::iterate(int n)
{
  int eflag,vflag,nflag,ntimestep;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;

    // initial time integration

    modify->initial_integrate();

    // regular communication vs neighbor list rebuild

    nflag = neighbor->decide();

    if (nflag == 0) {
      timer->stamp();
      comm->communicate();
      timer->stamp(TIME_COMM);
    } else {
      if (modify->n_pre_exchange) modify->pre_exchange();
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
      if (modify->n_pre_neighbor) modify->pre_neighbor();
      neighbor->build();
      timer->stamp(TIME_NEIGHBOR);
    }

    // eflag/vflag = 0/1/2 for energy/virial computation

    if (ntimestep == output->next_thermo) eflag = 1;
    else eflag = 0;

    if (ntimestep == output->next_thermo || ntimestep == next_virial) {
      vflag = virial_style;
      if (virial_every) next_virial++;
      else next_virial = fix_virial(ntimestep);
    } else vflag = 0;

    // force computations

    force_clear(vflag);

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

    if (modify->n_post_force) modify->post_force(vflag);
    modify->final_integrate();
    if (modify->n_end_of_step) modify->end_of_step();

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
}

/* ----------------------------------------------------------------------
   return next timestep virial should be computed
   based on one or more fixes that need virial computed periodically
------------------------------------------------------------------------- */

int Verlet::fix_virial(int ntimestep)
{
  int next;
  for (int ivirial = 0; ivirial < nfix_virial; ivirial++) {
    if (ntimestep == next_fix_virial[ivirial])
      next_fix_virial[ivirial] += fix_virial_every[ivirial];
    if (ivirial) next = MIN(next,next_fix_virial[ivirial]);
    else next = next_fix_virial[0];
  }
  return next;
}
