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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "hyper.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "integrate.h"
#include "min.h"
#include "force.h"
#include "neighbor.h"
#include "modify.h"
#include "compute_event_displace.h"
#include "fix_hyper.h"
#include "fix_event_hyper.h"
#include "output.h"
#include "dump.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOHYPER,GLOBAL,LOCAL};

/* ---------------------------------------------------------------------- */

Hyper::Hyper(LAMMPS *lmp) :
  Pointers(lmp), dumplist(NULL)
{}

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

  int nsteps = force->inumeric(FLERR,arg[0]);
  t_event = force->inumeric(FLERR,arg[1]);

  char *id_fix = new char[strlen(arg[2])+1];
  strcpy(id_fix,arg[2]);

  char *id_compute = new char[strlen(arg[3])+1];
  strcpy(id_compute,arg[3]);

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
    fix_hyper = (FixHyper *) modify->fix[ifix];
    int dim;
    int *hyperflag = (int *) fix_hyper->extract("hyperflag",dim);
    if (hyperflag == NULL || *hyperflag == 0)
      error->all(FLERR,"Hyper fix is not a valid hyperdynamics fix");
    if (*hyperflag == 1) hyperstyle = GLOBAL;
    if (*hyperflag == 2) hyperstyle = LOCAL;
    hyperenable = 1;
  }

  // create FixEventHyper class to store event and pre-quench states

  char **args = new char*[3];
  args[0] = (char *) "hyper_event";
  args[1] = (char *) "all";
  args[2] = (char *) "EVENT/HYPER";
  modify->add_fix(3,args);
  fix_event = (FixEventHyper *) modify->fix[modify->nfix-1];
  delete [] args;

  // create Finish for timing output

  finish = new Finish(lmp);

  // assign FixEventHyper to event-detection compute
  // necessary so it will know atom coords at last event

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID for hyper");
  compute_event = (ComputeEventDisplace *) modify->compute[icompute];
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
  update->endstep = update->laststep = update->firststep + nsteps;
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

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->time_depend)
      error->all(FLERR,"Cannot use hyper with a time-dependent fix defined");

  for (int i = 0; i < domain->nregion; i++)
    if (domain->regions[i]->dynamic_check())
      error->all(FLERR,"Cannot use hyper with a time-dependent region defined");

  // perform hyperdynamics simulation

  timer->init();
  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);

  nbuild = ndanger = 0;
  time_dynamics = time_quench = 0.0;

  if (hyperenable) fix_hyper->init_hyper();

  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);

  // perform initial minimization and bond list creation

  int nevent = 0;
  int nevent_atoms = 0;

  fix_event->store_state_quench();
  quench(1);
  if (dumpflag)
    for (int idump = 0; idump < ndump; idump++)
      output->dump[dumplist[idump]]->write();
  fix_event->store_event();
  if (hyperenable) fix_hyper->build_bond_list(0);
  fix_event->restore_state_quench();

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

      if (dumpflag)
        for (int idump = 0; idump < ndump; idump++)
          output->dump[dumplist[idump]]->write();
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

  if (me == 0) {
    if (screen) fprintf(screen,"Final hyper stats ...\n\n");
    if (logfile) fprintf(logfile,"Final hyper stats ...\n\n");
  }

  // subset of quantities also available in fix hyper output

  int nevent_running = 0;
  int nevent_atoms_running = 0;
  double t_hyper = 0.0;
  double avebonds = 0.0;
  double maxdrift = 0.0;
  double maxbondlen = 0.0;
  double fraczero = 1.0;

  double nnewbond,avenboost,aveboostcoeff,maxboostcoeff,minboostcoeff;
  double maxbondperatom,neighbondperbond,aveboostnow;
  double tbondbuild,rmaxever,rmaxeverbig,allghost_toofar;
  double lostbond,lostbondcoeff,biasoverlap,nonmatchbiascoeff;

  if (hyperenable) {
    t_hyper = fix_hyper->query(1);
    nevent_running = fix_hyper->query(2);
    nevent_atoms_running = fix_hyper->query(3);
    avebonds = fix_hyper->query(4);
    maxdrift = fix_hyper->query(5);
    maxbondlen = fix_hyper->query(6);
    fraczero = fix_hyper->query(7);

    if (hyperstyle == LOCAL) {
      nnewbond = fix_hyper->query(8);
      maxbondperatom = fix_hyper->query(9);
      avenboost = fix_hyper->query(10);
      aveboostcoeff = fix_hyper->query(11);
      maxboostcoeff = fix_hyper->query(12);
      minboostcoeff = fix_hyper->query(13);
      neighbondperbond = fix_hyper->query(14);
      aveboostnow = fix_hyper->query(15);
      tbondbuild = fix_hyper->query(16);
      rmaxever = fix_hyper->query(17);
      rmaxeverbig = fix_hyper->query(18);
      allghost_toofar = fix_hyper->query(19);
      lostbond = fix_hyper->query(20);
      lostbondcoeff = fix_hyper->query(21);
      biasoverlap = fix_hyper->query(22);
      nonmatchbiascoeff = fix_hyper->query(23);
    }
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Cummulative quantities for fix hyper:\n");
      fprintf(screen,"  hyper time = %g\n",t_hyper);
      fprintf(screen,"  time boost factor = %g\n",t_hyper/(nsteps*update->dt));
      fprintf(screen,"  event timesteps = %d\n",nevent_running);
      fprintf(screen,"  # of atoms in events = %d\n",nevent_atoms_running);
      fprintf(screen,"Quantities for this hyper run:\n");
      fprintf(screen,"  event timesteps = %d\n",nevent);
      fprintf(screen,"  # of atoms in events = %d\n",nevent_atoms);
      fprintf(screen,"  max length of any bond = %g\n",maxbondlen);
      fprintf(screen,"  max drift distance of any atom = %g\n",maxdrift);
      fprintf(screen,"  fraction of steps & bonds with zero bias = %g\n",
              fraczero);
      fprintf(screen,"Current quantities:\n");
      fprintf(screen,"  ave bonds/atom = %g\n",avebonds);

      if (hyperstyle == LOCAL) {
        fprintf(screen,"Cummulative quantities specific to fix hyper/local:\n");
        fprintf(screen,"  # of new bonds formed = %g\n",nnewbond);
        fprintf(screen,"  max bonds/atom = %g\n",maxbondperatom);
        fprintf(screen,"Quantities for this hyper run specific to "
                "fix hyper/local:\n");
        fprintf(screen,"  ave boosted bonds/step = %g\n",avenboost);
        fprintf(screen,"  ave boost coeff of all bonds = %g\n",aveboostcoeff);
        fprintf(screen,"  max boost coeff of any bond = %g\n",maxboostcoeff);
        fprintf(screen,"  min boost coeff of any bond = %g\n",minboostcoeff);
        fprintf(screen,"  max dist from my box of any "
                "non-maxstrain bond ghost atom = %g\n",rmaxever);
        fprintf(screen,"  max dist from my box of any bond ghost atom = %g\n",
                rmaxeverbig);
        fprintf(screen,"  count of bond ghost neighbors "
                "not found on reneighbor steps = %g\n",allghost_toofar);
        fprintf(screen,"  lost bond partners = %g\n",lostbond);
        fprintf(screen,"  ave bias coeff for lost bond partners = %g\n",
                lostbondcoeff);
        fprintf(screen,"  bias overlaps = %g\n",biasoverlap);
        fprintf(screen,"  non-matching bias coeffs = %g\n",nonmatchbiascoeff);
        fprintf(screen,"  CPU time for bond builds = %g\n",tbondbuild);
        fprintf(screen,"Current quantities specific to fix hyper/local:\n");
        fprintf(screen,"  neighbor bonds/bond = %g\n",neighbondperbond);
        fprintf(screen,"  ave boost coeff for all bonds = %g\n",aveboostnow);
      }
      fprintf(screen,"\n");
    }

    if (logfile) {
      fprintf(logfile,"Cummulative quantities for fix hyper:\n");
      fprintf(logfile,"  hyper time = %g\n",t_hyper);
      fprintf(logfile,"  event timesteps = %d\n",nevent_running);
      fprintf(logfile,"  # of atoms in events = %d\n",nevent_atoms_running);
      fprintf(logfile,"Quantities for this hyper run:\n");
      fprintf(logfile,"  event timesteps = %d\n",nevent);
      fprintf(logfile,"  # of atoms in events = %d\n",nevent_atoms);
      fprintf(logfile,"  max length of any bond = %g\n",maxbondlen);
      fprintf(logfile,"  max drift distance of any atom = %g\n",maxdrift);
      fprintf(logfile,"  fraction of steps & bonds with zero bias = %g\n",
              fraczero);
      fprintf(logfile,"Current quantities:\n");
      fprintf(logfile,"  ave bonds/atom = %g\n",avebonds);

      if (hyperstyle == LOCAL) {
        fprintf(logfile,"Cummulative quantities specific tofix hyper/local:\n");
        fprintf(logfile,"  # of new bonds formed = %g\n",nnewbond);
        fprintf(logfile,"  max bonds/atom = %g\n",maxbondperatom);
        fprintf(logfile,"Quantities for this hyper run specific to "
                "fix hyper/local:\n");
        fprintf(logfile,"  ave boosted bonds/step = %g\n",avenboost);
        fprintf(logfile,"  ave boost coeff of all bonds = %g\n",aveboostcoeff);
        fprintf(logfile,"  max boost coeff of any bond = %g\n",maxboostcoeff);
        fprintf(logfile,"  min boost coeff of any bond = %g\n",minboostcoeff);
        fprintf(logfile,"  max dist from my box of any "
                "non-maxstrain bond ghost atom = %g\n",rmaxever);
        fprintf(logfile,"  max dist from my box of any bond ghost atom = %g\n",
                rmaxeverbig);
        fprintf(logfile,"  count of ghost bond neighbors "
                "not found on reneighbor steps = %g\n",allghost_toofar);
        fprintf(logfile,"  lost bond partners = %g\n",lostbond);
        fprintf(logfile,"  ave bias coeff for lost bond partners = %g\n",
                lostbondcoeff);
        fprintf(logfile,"  bias overlaps = %g\n",biasoverlap);
        fprintf(logfile,"  non-matching bias coeffs = %g\n",nonmatchbiascoeff);
        fprintf(logfile,"  CPU time for bond builds = %g\n",tbondbuild);
        fprintf(logfile,"Current quantities specific to fix hyper/local:\n");
        fprintf(logfile,"  neighbor bonds/bond = %g\n",neighbondperbond);
        fprintf(logfile,"  ave boost coeff for all bonds = %g\n",aveboostnow);
      }
      fprintf(logfile,"\n");
    }
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
  memory->destroy(dumplist);
  delete finish;
  modify->delete_fix("hyper_event");

  compute_event->reset_extra_compute_fix(NULL);
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
  update->endstep = update->laststep = update->firststep + maxiter;
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
  ndump = 0;
  dumplist = NULL;
  rebond = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"min") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal hyper command");
      etol = force->numeric(FLERR,arg[iarg+1]);
      ftol = force->numeric(FLERR,arg[iarg+2]);
      maxiter = force->inumeric(FLERR,arg[iarg+3]);
      maxeval = force->inumeric(FLERR,arg[iarg+4]);
      if (maxiter < 0) error->all(FLERR,"Illegal hyper command");
      iarg += 5;

    } else if (strcmp(arg[iarg],"dump") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal hyper command");
      dumpflag = 1;
      int idump = output->find_dump(arg[iarg+1]);
      if (idump < 0)
        error->all(FLERR,"Dump ID in hyper command does not exist");
      memory->grow(dumplist,ndump+1,"hyper:dumplist");
      dumplist[ndump++] = idump;
      iarg += 2;

    } else if (strcmp(arg[iarg],"rebond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal hyper command");
      rebond = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal hyper command");
  }
}
