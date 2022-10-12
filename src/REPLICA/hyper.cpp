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

#include "hyper.h"

#include "compute_event_displace.h"
#include "domain.h"
#include "dump.h"
#include "error.h"
#include "finish.h"
#include "fix_event_hyper.h"
#include "fix_hyper.h"
#include "integrate.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "region.h"
#include "timer.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{NOHYPER,GLOBAL,LOCAL};

/* ---------------------------------------------------------------------- */

Hyper::Hyper(LAMMPS *_lmp) : Command(_lmp) {}

/* ----------------------------------------------------------------------
   perform hyperdynamics simulation
------------------------------------------------------------------------- */

void Hyper::command(int narg, char **arg)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // error checks

  if (domain->box_exist == 0)
    error->all(FLERR,"Hyper command before simulation box is defined");

  if (narg < 4) error->all(FLERR,"Illegal hyper command");

  int nsteps = utils::inumeric(FLERR,arg[0],false,lmp);
  t_event = utils::inumeric(FLERR,arg[1],false,lmp);

  auto id_fix = utils::strdup(arg[2]);
  auto id_compute = utils::strdup(arg[3]);

  options(narg-4,&arg[4]);

  // total # of timesteps must be multiple of t_event

  if (t_event <= 0)
    error->all(FLERR,"Invalid t_event in hyper command");
  if (nsteps % t_event)
    error->all(FLERR,"Hyper nsteps must be multiple of t_event");
  if (rebond < 0)
    error->all(FLERR,"Invalid rebond in hyper command");
  if (rebond && rebond % t_event)
    error->all(FLERR,"Hyper rebond must be multiple of t_event");

  // FixHyper class performs global or local hyperdynamics

  int hyperenable,hyperstyle;

  if (strcmp(id_fix,"NULL") == 0) {
    hyperenable = 0;
    hyperstyle = NOHYPER;
  } else {
    int ifix = modify->find_fix(id_fix);
    if (ifix < 0) error->all(FLERR,"Could not find fix ID for hyper");
    fix_hyper = dynamic_cast<FixHyper *>(modify->fix[ifix]);
    int dim;
    int *hyperflag = (int *) fix_hyper->extract("hyperflag",dim);
    if (hyperflag == nullptr || *hyperflag == 0)
      error->all(FLERR,"Hyper fix is not a valid hyperdynamics fix");
    if (*hyperflag == 1) hyperstyle = GLOBAL;
    if (*hyperflag == 2) hyperstyle = LOCAL;
    hyperenable = 1;
  }

  // create FixEventHyper class to store event and pre-quench states

  fix_event = dynamic_cast<FixEventHyper *>(modify->add_fix("hyper_event all EVENT/HYPER"));

  // create Finish for timing output

  finish = new Finish(lmp);

  // assign FixEventHyper to event-detection compute
  // necessary so it will know atom coords at last event

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID for hyper");
  compute_event = dynamic_cast<ComputeEventDisplace *>(modify->compute[icompute]);
  compute_event->reset_extra_compute_fix("hyper_event");

  // reset reneighboring criteria since will perform minimizations

  neigh_every = neighbor->every;
  neigh_delay = neighbor->delay;
  neigh_dist_check = neighbor->dist_check;

  if (neigh_every != 1 || neigh_delay != 0 || neigh_dist_check != 1) {
    if (me == 0)
      error->warning(FLERR,"Resetting reneighboring criteria during hyper");
  }

  neighbor->every = 1;
  neighbor->delay = 0;
  neighbor->dist_check = 1;

  // initialize hyper as if one long dynamics run

  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->beginstep + nsteps;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // init minimizer settings and minimizer itself

  update->etol = etol;
  update->ftol = ftol;
  update->max_eval = maxeval;

  // cannot use hyper with a changing box
  // removing this restriction would require saving/restoring box params

  if (domain->box_change)
    error->all(FLERR,"Cannot use hyper with a changing box");

  // cannot use hyper with time-dependent fixes or regions

  for (auto &ifix : modify->get_fix_list())
    if (ifix->time_depend) error->all(FLERR,"Cannot use hyper with a time-dependent fix defined");

  for (auto &reg : domain->get_region_list())
    if (reg->dynamic_check())
      error->all(FLERR,"Cannot use hyper with a time-dependent region defined");

  // perform hyperdynamics simulation

  timer->init();
  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);

  nbuild = ndanger = 0;
  time_dynamics = time_quench = 0.0;

  if (hyperenable) fix_hyper->init_hyper();

  // perform initial minimization and bond list creation

  int nevent = 0;
  int nevent_atoms = 0;

  fix_event->store_state_quench();
  quench(1);
  if (dumpflag) for (auto &idump : dumplist) idump->write();
  fix_event->store_event();
  if (hyperenable) fix_hyper->build_bond_list(0);
  fix_event->restore_state_quench();

  // reset stats and timers to skip HD setup

  nbuild = ndanger = 0;
  time_dynamics = time_quench = 0.0;

  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);

  // main loop: dynamics, store state, quench, check event, restore state

  int ecount;
  int istep = 0;

  while (istep < nsteps) {
    dynamics(t_event,time_dynamics);
    fix_event->store_state_quench();
    quench(0);

    ecount = compute_event->all_events();

    if (ecount) {
      nevent++;
      nevent_atoms += ecount;

      if (dumpflag) for (auto &idump : dumplist) idump->write();
      fix_event->store_event();
      if (hyperenable) fix_hyper->build_bond_list(ecount);

    } else if (rebond && update->ntimestep % rebond == 0) {
      fix_event->store_event();
      if (hyperenable) fix_hyper->build_bond_list(ecount);
    }

    fix_event->restore_state_quench();
    istep = update->ntimestep - update->beginstep;
  }

  nsteps = update->ntimestep - update->beginstep;

  // set total timers and counters so Finish() will process them

  timer->set_wall(Timer::TOTAL,time_start);
  timer->barrier_stop();

  timer->set_wall(Timer::DYNAMICS,time_dynamics);
  timer->set_wall(Timer::QUENCH,time_quench);

  neighbor->ncalls = nbuild;
  neighbor->ndanger = ndanger;

  update->nsteps = nsteps;

  if (me == 0) utils::logmesg(lmp,"Final hyper stats ...\n\n");

  // subset of quantities also available in fix hyper output
  // set t_hyper to no-boost value when hyperenable is not set

  int nevent_running = 0;
  int nevent_atoms_running = 0;
  double t_hyper = update->dt * (update->endstep - update->beginstep);
  double avebonds = 0.0;
  double maxdrift = 0.0;
  double maxbondlen = 0.0;
  double fraczero = 1.0;
  double fracneg = 1.0;

  double nnewbond,aveboost,avenbias,avebiascoeff,minbiascoeff,maxbiascoeff;
  double maxbondperatom,neighbondperbond,avebiasnow;
  double tbondbuild,rmaxever,rmaxeverbig,allghost_toofar;
  double biasoverlap;

  if (hyperenable) {
    t_hyper = fix_hyper->query(1);
    nevent_running = fix_hyper->query(2);
    nevent_atoms_running = fix_hyper->query(3);
    avebonds = fix_hyper->query(4);
    maxdrift = fix_hyper->query(5);
    maxbondlen = fix_hyper->query(6);
    fraczero = fix_hyper->query(7);
    fracneg = fix_hyper->query(8);

    if (hyperstyle == LOCAL) {
      nnewbond = fix_hyper->query(9);
      maxbondperatom = fix_hyper->query(10);
      aveboost = fix_hyper->query(11);
      avenbias = fix_hyper->query(12);
      avebiascoeff = fix_hyper->query(13);
      minbiascoeff = fix_hyper->query(14);
      maxbiascoeff = fix_hyper->query(15);
      neighbondperbond = fix_hyper->query(16);
      avebiasnow = fix_hyper->query(17);
      tbondbuild = fix_hyper->query(18);
      rmaxever = fix_hyper->query(19);
      rmaxeverbig = fix_hyper->query(20);
      allghost_toofar = fix_hyper->query(21);
      biasoverlap = fix_hyper->query(22);
    }
  }

  if (me == 0) {
    std::string mesg = "Cummulative quantities for fix hyper:\n";
    mesg += fmt::format("  hyper time = {}\n",t_hyper);
    if (hyperenable)
      mesg += fmt::format("  time boost factor = {}\n", t_hyper /
                          ((update->ntimestep -fix_hyper->ntimestep_initial)*update->dt));
    else mesg += "  time boost factor = 1\n";
    mesg += fmt::format("  event timesteps = {}\n",nevent_running);
    mesg += fmt::format("  # of atoms in events = {}\n",nevent_atoms_running);
    mesg += "Quantities for this hyper run:\n";
    mesg += fmt::format("  event timesteps = {}\n",nevent);
    mesg += fmt::format("  # of atoms in events = {}\n",nevent_atoms);
    mesg += fmt::format("  max length of any bond = {}\n",maxbondlen);
    mesg += fmt::format("  max drift distance of any atom = {}\n",maxdrift);
    mesg += fmt::format("  fraction of biased bonds with zero bias = {}\n",fraczero);
    mesg += fmt::format("  fraction of biased bonds with negative strain = {}\n",fracneg);
    mesg += "Current quantities:\n";
    mesg += fmt::format("  ave bonds/atom = {}\n",avebonds);

    if (hyperstyle == LOCAL) {
      mesg += "Cummulative quantities specific to fix hyper/local:\n";
      mesg += fmt::format("  # of new bonds formed = {}\n",nnewbond);
      mesg += fmt::format("  max bonds/atom = {}\n",maxbondperatom);
      mesg += "Quantities for this hyper run specific to fix hyper/local:\n";
      mesg += fmt::format("  ave boost for all bonds/step = {}\n",aveboost);
      mesg += fmt::format("  ave biased bonds/step = {}\n",avenbias);
      mesg += fmt::format("  ave bias coeff of all bonds = {}\n",avebiascoeff);
      mesg += fmt::format("  min bias coeff of any bond = {}\n",minbiascoeff);
      mesg += fmt::format("  max bias coeff of any bond = {}\n",maxbiascoeff);
      mesg += fmt::format("  max dist from my subbox of any "
                          "non-maxstrain bond ghost atom = {}\n",rmaxever);
      mesg += fmt::format("  max dist from my box of any bond ghost atom = {}\n",
                          rmaxeverbig);
      mesg += fmt::format("  count of bond ghost neighbors "
                          "not found on reneighbor steps = {}\n",allghost_toofar);
      mesg += fmt::format("  bias overlaps = {}\n",biasoverlap);
      mesg += fmt::format("  CPU time for bond builds = {}\n",tbondbuild);
      mesg += "Current quantities specific to fix hyper/local:\n";
      mesg += fmt::format("  neighbor bonds/bond = {}\n",neighbondperbond);
      mesg += fmt::format("  ave boost coeff for all bonds = {}\n",avebiasnow);
    }
    utils::logmesg(lmp, mesg);
  }

  // timing stats

  finish->end(4);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;

  // reset reneighboring criteria

  neighbor->every = neigh_every;
  neighbor->delay = neigh_delay;
  neighbor->dist_check = neigh_dist_check;

  delete [] id_fix;
  delete [] id_compute;
  delete finish;
  modify->delete_fix("hyper_event");

  compute_event->reset_extra_compute_fix(nullptr);
}

/* ----------------------------------------------------------------------
   short dynamics run
------------------------------------------------------------------------- */

void Hyper::dynamics(int nsteps, double & /* time_category */)
{
  update->whichflag = 1;
  update->nsteps = nsteps;

  // full init works
  // need to try partial init or setup

  lmp->init();
  update->integrate->setup(0);

  // this may be needed if don't do full init
  //modify->addstep_compute_all(update->ntimestep);
  bigint ncalls = neighbor->ncalls;

  timer->barrier_start();
  update->integrate->run(nsteps);
  timer->barrier_stop();
  time_dynamics += timer->get_wall(Timer::TOTAL);

  nbuild += neighbor->ncalls - ncalls;
  ndanger += neighbor->ndanger;

  update->integrate->cleanup();
  finish->end(0);
}

/* ----------------------------------------------------------------------
   quench minimization
   flag = 1 to trigger output of memory in setup() call
------------------------------------------------------------------------- */

void Hyper::quench(int flag)
{
  bigint ntimestep_hold = update->ntimestep;
  bigint endstep_hold = update->endstep;

  // need to change whichflag so that minimize->setup() calling
  // modify->setup() will call fix->min_setup()

  update->whichflag = 2;
  update->nsteps = maxiter;
  update->endstep = update->laststep = update->ntimestep + maxiter;
  if (update->laststep < 0)
    error->all(FLERR,"Too many iterations");
  update->restrict_output = 1;

  // full init works

  lmp->init();
  update->minimize->setup(flag);

  // partial init does not work

  //modify->addstep_compute_all(update->ntimestep);
  //update->minimize->setup_minimal(1);

  timer->barrier_start();
  update->minimize->run(maxiter);
  timer->barrier_stop();
  time_quench += timer->get_wall(Timer::TOTAL);

  update->minimize->cleanup();
  finish->end(0);

  // reset timestep as if quench did not occur
  // clear timestep storage from computes, since now invalid

  update->restrict_output = 0;
  update->ntimestep = ntimestep_hold;
  update->endstep = update->laststep = endstep_hold;
  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of hyper input line
------------------------------------------------------------------------- */

void Hyper::options(int narg, char **arg)
{
  // set defaults

  etol = 1.0e-4;
  ftol = 1.0e-4;
  maxiter = 40;
  maxeval = 50;
  dumpflag = 0;
  rebond = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"min") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal hyper command");
      etol = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ftol = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      maxiter = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      maxeval = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      if (maxiter < 0) error->all(FLERR,"Illegal hyper command");
      iarg += 5;

    } else if (strcmp(arg[iarg],"dump") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal hyper command");
      dumpflag = 1;
      auto idump = output->get_dump_by_id(arg[iarg+1]);
      if (!idump) error->all(FLERR,"Dump ID {} in hyper command does not exist", arg[iarg+1]);
      dumplist.emplace_back(idump);
      iarg += 2;

    } else if (strcmp(arg[iarg],"rebond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal hyper command");
      rebond = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else error->all(FLERR,"Illegal hyper command");
  }
}
