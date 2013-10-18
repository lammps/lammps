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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "prd.h"
#include "universe.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "velocity.h"
#include "integrate.h"
#include "min.h"
#include "neighbor.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_event_prd.h"
#include "force.h"
#include "pair.h"
#include "random_park.h"
#include "random_mars.h"
#include "output.h"
#include "dump.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PRD::PRD(LAMMPS *lmp) : Pointers(lmp) {}

/* ----------------------------------------------------------------------
   perform PRD simulation on one or more replicas
------------------------------------------------------------------------- */

void PRD::command(int narg, char **arg)
{
  int flag,ireplica;

  // error checks

  if (domain->box_exist == 0)
    error->all(FLERR,"PRD command before simulation box is defined");
  if (universe->nworlds != universe->nprocs &&
      atom->map_style == 0)
    error->all(FLERR,"Cannot use PRD with multi-processor replicas "
               "unless atom map exists");
  if (universe->nworlds == 1 && comm->me == 0)
    error->warning(FLERR,"Running PRD with only one replica");

  if (narg < 7) error->universe_all(FLERR,"Illegal prd command");

  nsteps = force->inumeric(FLERR,arg[0]);
  t_event = force->inumeric(FLERR,arg[1]);
  n_dephase = force->inumeric(FLERR,arg[2]);
  t_dephase = force->inumeric(FLERR,arg[3]);
  t_corr = force->inumeric(FLERR,arg[4]);

  char *id_compute = new char[strlen(arg[5])+1];
  strcpy(id_compute,arg[5]);
  int seed = force->inumeric(FLERR,arg[6]);

  options(narg-7,&arg[7]);

  // total # of timesteps must be multiple of t_event

  if (t_event <= 0) error->universe_all(FLERR,"Invalid t_event in prd command");
  if (nsteps % t_event)
    error->universe_all(FLERR,"PRD nsteps must be multiple of t_event");
  if (t_corr % t_event)
    error->universe_all(FLERR,"PRD t_corr must be multiple of t_event");

  // local storage

  int me_universe = universe->me;
  int nprocs_universe = universe->nprocs;
  int nreplica = universe->nworlds;
  int iworld = universe->iworld;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // comm_replica = communicator between all proc 0s across replicas

  int color = me;
  MPI_Comm_split(universe->uworld,color,0,&comm_replica);

  // equal_size_replicas = 1 if all replicas have same # of procs
  // no longer used

  //flag = 0;
  //if (nreplica*nprocs == nprocs_universe) flag = 1;
  //MPI_Allreduce(&flag,&equal_size_replicas,1,MPI_INT,MPI_MIN,
  //              universe->uworld);

  // workspace for inter-replica communication via gathers

  natoms = atom->natoms;

  displacements = NULL;
  tagall = NULL;
  xall = NULL;
  imageall = NULL;

  if (nreplica != nprocs_universe) {
    displacements = new int[nprocs];
    memory->create(tagall,natoms,"prd:tagall");
    memory->create(xall,natoms,3,"prd:xall");
    memory->create(imageall,natoms,"prd:imageall");
  }

  // random_select = same RNG for each replica for multiple event selection
  // random_dephase = unique RNG for each replica for dephasing

  random_select = new RanPark(lmp,seed);
  random_dephase = new RanMars(lmp,seed+iworld);

  // create ComputeTemp class to monitor temperature

  char **args = new char*[3];
  args[0] = (char *) "prd_temp";
  args[1] = (char *) "all";
  args[2] = (char *) "temp";
  modify->add_compute(3,args);
  temperature = modify->compute[modify->ncompute-1];

  // create Velocity class for velocity creation in dephasing
  // pass it temperature compute, loop_setting, dist_setting settings

  atom->check_mass();
  velocity = new Velocity(lmp);
  velocity->init_external("all");

  args[0] = (char *) "temp";
  args[1] = (char *) "prd_temp";
  velocity->options(2,args);
  args[0] = (char *) "loop";
  args[1] = (char *) loop_setting;
  if (loop_setting) velocity->options(2,args);
  args[0] = (char *) "dist";
  args[1] = (char *) dist_setting;
  if (dist_setting) velocity->options(2,args);

  // create FixEventPRD class to store event and pre-quench states

  args[0] = (char *) "prd_event";
  args[1] = (char *) "all";
  args[2] = (char *) "EVENT/PRD";
  modify->add_fix(3,args);
  fix_event = (FixEventPRD *) modify->fix[modify->nfix-1];

  // create Finish for timing output

  finish = new Finish(lmp);

  // string clean-up

  delete [] args;
  delete [] loop_setting;
  delete [] dist_setting;

  // assign FixEventPRD to event-detection compute
  // necessary so it will know atom coords at last event

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID for PRD");
  compute_event = modify->compute[icompute];
  compute_event->reset_extra_compute_fix("prd_event");

  // reset reneighboring criteria since will perform minimizations

  neigh_every = neighbor->every;
  neigh_delay = neighbor->delay;
  neigh_dist_check = neighbor->dist_check;

  if (neigh_every != 1 || neigh_delay != 0 || neigh_dist_check != 1) {
    if (me == 0)
      error->warning(FLERR,"Resetting reneighboring criteria during PRD");
  }

  neighbor->every = 1;
  neighbor->delay = 0;
  neighbor->dist_check = 1;

  // initialize PRD as if one long dynamics run

  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  update->restrict_output = 1;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // init minimizer settings and minimizer itself

  update->etol = etol;
  update->ftol = ftol;
  update->max_eval = maxeval;

  update->minimize->init();

  // cannot use PRD with a changing box

  if (domain->box_change)
    error->all(FLERR,"Cannot use PRD with a changing box");

  // cannot use PRD with time-dependent fixes or regions or atom sorting

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->time_depend)
      error->all(FLERR,"Cannot use PRD with a time-dependent fix defined");

  for (int i = 0; i < domain->nregion; i++)
    if (domain->regions[i]->dynamic_check())
      error->all(FLERR,"Cannot use PRD with a time-dependent region defined");

  if (atom->sortfreq > 0)
    error->all(FLERR,"Cannot use PRD with atom_modify sort enabled");

  // perform PRD simulation

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up PRD ...\n");

  if (me_universe == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,"Step CPU Clock Event "
              "Correlated Coincident Replica\n");
    if (universe->ulogfile)
      fprintf(universe->ulogfile,"Step CPU Clock Event "
              "Correlated Coincident Replica\n");
  }

  // store hot state and quenched event for replica 0
  // use share_event() to copy that info to all replicas
  // this insures all start from same place

  // need this line if quench() does only setup_minimal()
  // update->minimize->setup();

  fix_event->store_state();
  quench();
  ncoincident = 0;
  share_event(0,0);

  timer->init();
  timer->barrier_start(TIME_LOOP);
  time_start = timer->array[TIME_LOOP];

  log_event();

  // do full init/setup since are starting all replicas after event
  // replica 0 bcasts temp to all replicas if temp_dephase is not set

  update->whichflag = 1;
  lmp->init();
  update->integrate->setup();

  if (temp_flag == 0) {
    if (universe->iworld == 0) temp_dephase = temperature->compute_scalar();
    MPI_Bcast(&temp_dephase,1,MPI_DOUBLE,universe->root_proc[0],
              universe->uworld);
  }

  // main loop: look for events until out of time
  // (1) dephase independently on each proc after event
  // (2) loop: dynamics, store state, quench, check event, restore state
  // (3) share and record event

  nbuild = ndanger = 0;
  time_dephase = time_dynamics = time_quench = time_comm = time_output = 0.0;

  timer->barrier_start(TIME_LOOP);
  time_start = timer->array[TIME_LOOP];

  while (update->ntimestep < update->endstep) {
    dephase();

    ireplica = -1;
    while (update->ntimestep < update->endstep) {
      dynamics();
      fix_event->store_state();
      quench();
      ireplica = check_event();
      if (ireplica >= 0) break;
      fix_event->restore_state();
    }
    if (ireplica < 0) break;

    // potentially more efficient for correlated events if don't
    // share until correlated check has completed
    // this will complicate the dump (always on replica 0)

    share_event(ireplica,1);
    log_event();

    int restart_flag = 0;
    if (output->restart_flag && universe->iworld == 0) {
      if (output->restart_every_single &&
          fix_event->event_number % output->restart_every_single == 0)
        restart_flag = 1;
      if (output->restart_every_double &&
          fix_event->event_number % output->restart_every_double == 0)
        restart_flag = 1;
    }

    // correlated event loop
    // other procs could be dephasing during this time

    int corr_endstep = update->ntimestep + t_corr;
    while (update->ntimestep < corr_endstep) {
      if (update->ntimestep == update->endstep) {
        restart_flag = 0;
        break;
      }
      dynamics();
      fix_event->store_state();
      quench();
      int corr_event_check = check_event(ireplica);
      if (corr_event_check >= 0) {
        share_event(ireplica,2);
        log_event();
        corr_endstep = update->ntimestep + t_corr;
      } else fix_event->restore_state();
    }

    // full init/setup since are starting all replicas after event
    // event replica bcasts temp to all replicas if temp_dephase is not set

    update->whichflag = 1;
    lmp->init();
    update->integrate->setup();

    timer->barrier_start(TIME_LOOP);

    if (t_corr > 0) replicate(ireplica);
    if (temp_flag == 0) {
      if (ireplica == universe->iworld)
        temp_dephase = temperature->compute_scalar();
      MPI_Bcast(&temp_dephase,1,MPI_DOUBLE,universe->root_proc[ireplica],
                      universe->uworld);
    }

    timer->barrier_stop(TIME_LOOP);
    time_comm += timer->array[TIME_LOOP];

    // write restart file of hot coords

    if (restart_flag) {
      timer->barrier_start(TIME_LOOP);
      output->write_restart(update->ntimestep);
      timer->barrier_stop(TIME_LOOP);
      time_output += timer->array[TIME_LOOP];
    }
  }

  // set total timers and counters so Finish() will process them

  timer->array[TIME_LOOP] = time_start;
  timer->barrier_stop(TIME_LOOP);

  timer->array[TIME_PAIR] = time_dephase;
  timer->array[TIME_BOND] = time_dynamics;
  timer->array[TIME_KSPACE] = time_quench;
  timer->array[TIME_COMM] = time_comm;
  timer->array[TIME_OUTPUT] = time_output;

  neighbor->ncalls = nbuild;
  neighbor->ndanger = ndanger;

  if (me_universe == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,
              "Loop time of %g on %d procs for %d steps with " BIGINT_FORMAT
              " atoms\n",
              timer->array[TIME_LOOP],nprocs_universe,nsteps,atom->natoms);
    if (universe->ulogfile)
      fprintf(universe->ulogfile,
              "Loop time of %g on %d procs for %d steps with " BIGINT_FORMAT
              " atoms\n",
              timer->array[TIME_LOOP],nprocs_universe,nsteps,atom->natoms);
  }

  finish->end(2);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
  update->restrict_output = 0;

  // reset reneighboring criteria

  neighbor->every = neigh_every;
  neighbor->delay = neigh_delay;
  neighbor->dist_check = neigh_dist_check;

  // clean up

  delete [] displacements;
  memory->destroy(tagall);
  memory->destroy(xall);
  memory->destroy(imageall);

  delete [] id_compute;
  MPI_Comm_free(&comm_replica);
  delete random_select;
  delete random_dephase;
  delete velocity;
  delete finish;
  modify->delete_compute("prd_temp");
  modify->delete_fix("prd_event");

  compute_event->reset_extra_compute_fix(NULL);
}

/* ----------------------------------------------------------------------
   dephasing = one or more short runs with new random velocities
------------------------------------------------------------------------- */

void PRD::dephase()
{
  bigint ntimestep_hold = update->ntimestep;

  update->whichflag = 1;
  update->nsteps = n_dephase*t_dephase;

  timer->barrier_start(TIME_LOOP);

  for (int i = 0; i < n_dephase; i++) {
    int seed = static_cast<int> (random_dephase->uniform() * MAXSMALLINT);
    if (seed == 0) seed = 1;
    velocity->create(temp_dephase,seed);
    update->integrate->run(t_dephase);
    if (temp_flag == 0) temp_dephase = temperature->compute_scalar();
  }

  timer->barrier_stop(TIME_LOOP);
  time_dephase += timer->array[TIME_LOOP];

  update->integrate->cleanup();
  finish->end(0);

  // reset timestep as if dephase did not occur
  // clear timestep storage from computes, since now invalid

  update->ntimestep = ntimestep_hold;
  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();
}

/* ----------------------------------------------------------------------
   single short dynamics run
------------------------------------------------------------------------- */

void PRD::dynamics()
{
  update->whichflag = 1;
  update->nsteps = t_event;

  lmp->init();
  update->integrate->setup();
  // this may be needed if don't do full init
  //modify->addstep_compute_all(update->ntimestep);
  bigint ncalls = neighbor->ncalls;

  timer->barrier_start(TIME_LOOP);
  update->integrate->run(t_event);
  timer->barrier_stop(TIME_LOOP);
  time_dynamics += timer->array[TIME_LOOP];

  nbuild += neighbor->ncalls - ncalls;
  ndanger += neighbor->ndanger;

  update->integrate->cleanup();
  finish->end(0);
}

/* ----------------------------------------------------------------------
   quench minimization
------------------------------------------------------------------------- */

void PRD::quench()
{
  bigint ntimestep_hold = update->ntimestep;
  bigint endstep_hold = update->endstep;

  // need to change whichflag so that minimize->setup() calling
  // modify->setup() will call fix->min_setup()

  update->whichflag = 2;
  update->nsteps = maxiter;
  update->endstep = update->laststep = update->firststep + maxiter;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many iterations");

  // full init works

  lmp->init();
  update->minimize->setup();

  // partial init does not work

  //modify->addstep_compute_all(update->ntimestep);
  //update->minimize->setup_minimal(1);

  int ncalls = neighbor->ncalls;

  timer->barrier_start(TIME_LOOP);
  update->minimize->run(maxiter);
  timer->barrier_stop(TIME_LOOP);
  time_quench += timer->array[TIME_LOOP];

  if (neighbor->ncalls == ncalls) quench_reneighbor = 0;
  else quench_reneighbor = 1;

  update->minimize->cleanup();
  finish->end(0);

  // reset timestep as if dephase did not occur
  // clear timestep storage from computes, since now invalid

  update->ntimestep = ntimestep_hold;
  update->endstep = update->laststep = endstep_hold;
  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();
}

/* ----------------------------------------------------------------------
   check for an event in any replica
   if replica_num is non-negative only check for event on replica_num
   if multiple events, choose one at random
   return -1 if no event
   else return ireplica = world in which event occured
------------------------------------------------------------------------- */

int PRD::check_event(int replica_num)
{
  int worldflag,universeflag,scanflag,replicaflag,ireplica;

  worldflag = 0;
  if (compute_event->compute_scalar() > 0.0) worldflag = 1;
  if (replica_num >= 0 && replica_num != universe->iworld) worldflag = 0;

  timer->barrier_start(TIME_LOOP);

  if (me == 0) MPI_Allreduce(&worldflag,&universeflag,1,
                             MPI_INT,MPI_SUM,comm_replica);
  MPI_Bcast(&universeflag,1,MPI_INT,0,world);

  ncoincident = universeflag;

  if (!universeflag) ireplica = -1;
  else {

    // multiple events, choose one at random
    // iwhich = random # from 1 to N, N = # of events to choose from
    // scanflag = 1 to N on replicas with an event, 0 on non-event replicas
    // exit with worldflag = 1 on chosen replica, 0 on all others
    // note worldflag is already 0 on replicas that didn't perform event

    if (universeflag > 1) {
      int iwhich = static_cast<int>
        (universeflag*random_select->uniform()) + 1;

      if (me == 0)
        MPI_Scan(&worldflag,&scanflag,1,MPI_INT,MPI_SUM,comm_replica);
      MPI_Bcast(&scanflag,1,MPI_INT,0,world);

      if (scanflag != iwhich) worldflag = 0;
    }

    if (worldflag) replicaflag = universe->iworld;
    else replicaflag = 0;

    if (me == 0) MPI_Allreduce(&replicaflag,&ireplica,1,
                               MPI_INT,MPI_SUM,comm_replica);
    MPI_Bcast(&ireplica,1,MPI_INT,0,world);
  }

  timer->barrier_stop(TIME_LOOP);
  time_comm += timer->array[TIME_LOOP];

  return ireplica;
}

/* ----------------------------------------------------------------------
   share quenched and hot coords owned by ireplica with all replicas
   all replicas store event in fix_event
   replica 0 dumps event snapshot
   flag = 0 = called before PRD run
   flag = 1 = called during PRD run = not correlated event
   flag = 2 = called during PRD run = correlated event
------------------------------------------------------------------------- */

void PRD::share_event(int ireplica, int flag)
{
  timer->barrier_start(TIME_LOOP);

  // communicate quenched coords to all replicas and store as event
  // decrement event counter if flag = 0 since not really an event

  replicate(ireplica);
  timer->barrier_stop(TIME_LOOP);
  time_comm += timer->array[TIME_LOOP];

  // adjust time for last correlated event check (not on first event)

  int corr_adjust = t_corr;
  if (fix_event->event_number < 1 || flag == 2) corr_adjust = 0;

  // delta = time since last correlated event check

  int delta = update->ntimestep - fix_event->event_timestep - corr_adjust;

  // if this is a correlated event, time elapsed only on one partition

  if (flag != 2) delta *= universe->nworlds;
  delta += corr_adjust;

  // don't change the clock or timestep if this is a restart

  if (flag == 0 && fix_event->event_number != 0)
    fix_event->store_event_prd(fix_event->event_timestep,0);
  else {
    fix_event->store_event_prd(update->ntimestep,delta);
    fix_event->replica_number = ireplica;
    fix_event->correlated_event = 0;
    if (flag == 2) fix_event->correlated_event = 1;
    fix_event->ncoincident = ncoincident;
  }
  if (flag == 0) fix_event->event_number--;

  // dump snapshot of quenched coords, only on replica 0
  // must reneighbor and compute forces before dumping
  // since replica 0 possibly has new state from another replica
  // addstep_compute_all insures eng/virial are calculated if needed

  if (output->ndump && universe->iworld == 0) {
    timer->barrier_start(TIME_LOOP);
    modify->addstep_compute_all(update->ntimestep);
    update->integrate->setup_minimal(1);
    output->write_dump(update->ntimestep);
    timer->barrier_stop(TIME_LOOP);
    time_output += timer->array[TIME_LOOP];
  }

  // restore and communicate hot coords to all replicas

  fix_event->restore_state();
  timer->barrier_start(TIME_LOOP);
  replicate(ireplica);
  timer->barrier_stop(TIME_LOOP);
  time_comm += timer->array[TIME_LOOP];
}

/* ----------------------------------------------------------------------
   universe proc 0 prints event info
------------------------------------------------------------------------- */

void PRD::log_event()
{
  timer->array[TIME_LOOP] = time_start;
  if (universe->me == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,
              BIGINT_FORMAT " %.3f %d %d %d %d %d\n",
              fix_event->event_timestep,
              timer->elapsed(TIME_LOOP),
              fix_event->clock,
              fix_event->event_number,fix_event->correlated_event,
              fix_event->ncoincident,
              fix_event->replica_number);
    if (universe->ulogfile)
      fprintf(universe->ulogfile,
              BIGINT_FORMAT " %.3f %d %d %d %d %d\n",
              fix_event->event_timestep,
              timer->elapsed(TIME_LOOP),
              fix_event->clock,
              fix_event->event_number,fix_event->correlated_event,
              fix_event->ncoincident,
              fix_event->replica_number);
  }
}

/* ----------------------------------------------------------------------
  communicate atom coords and image flags in ireplica to all other replicas
  one proc per replica:
    direct overwrite via bcast
  else atoms could be stored in different order or on different procs:
    collect to root proc of event replica
    bcast to roots of other replicas
    bcast within each replica
    each proc extracts info for atoms it owns using atom IDs
------------------------------------------------------------------------- */

void PRD::replicate(int ireplica)
{
  int nreplica = universe->nworlds;
  int nprocs_universe = universe->nprocs;
  int i,m,flag,commflag;

  if (nreplica == nprocs_universe) {
    MPI_Bcast(atom->image,atom->nlocal,MPI_INT,ireplica,comm_replica);
    MPI_Bcast(atom->x[0],3*atom->nlocal,MPI_DOUBLE,ireplica,comm_replica);

  } else {
    int *counts = new int[nprocs];

    if (universe->iworld == ireplica) {
      MPI_Gather(&atom->nlocal,1,MPI_INT,counts,1,MPI_INT,0,world);
      displacements[0] = 0;
      for (i = 0; i < nprocs-1; i++)
        displacements[i+1] = displacements[i] + counts[i];
      MPI_Gatherv(atom->tag,atom->nlocal,MPI_INT,
                  tagall,counts,displacements,MPI_INT,0,world);
      MPI_Gatherv(atom->image,atom->nlocal,MPI_INT,
                        imageall,counts,displacements,MPI_INT,0,world);
      for (i = 0; i < nprocs; i++) counts[i] *= 3;
      for (i = 0; i < nprocs-1; i++)
        displacements[i+1] = displacements[i] + counts[i];
      MPI_Gatherv(atom->x[0],3*atom->nlocal,MPI_DOUBLE,
                        xall[0],counts,displacements,MPI_DOUBLE,0,world);
    }

    if (me == 0) {
      MPI_Bcast(tagall,natoms,MPI_INT,ireplica,comm_replica);
      MPI_Bcast(imageall,natoms,MPI_INT,ireplica,comm_replica);
      MPI_Bcast(xall[0],3*natoms,MPI_DOUBLE,ireplica,comm_replica);
    }

    MPI_Bcast(tagall,natoms,MPI_INT,0,world);
    MPI_Bcast(imageall,natoms,MPI_INT,0,world);
    MPI_Bcast(xall[0],3*natoms,MPI_DOUBLE,0,world);

    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (i = 0; i < natoms; i++) {
      m = atom->map(tagall[i]);
      if (m >= 0 && m < nlocal) {
        x[m][0] = xall[i][0];
        x[m][1] = xall[i][1];
        x[m][2] = xall[i][2];
        atom->image[m] = imageall[i];
      }
    }

    delete [] counts;
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of PRD input line
------------------------------------------------------------------------- */

void PRD::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal prd command");

  // set defaults

  etol = 0.1;
  ftol = 0.1;
  maxiter = 40;
  maxeval = 50;
  temp_flag = 0;

  char *str = (char *) "geom";
  int n = strlen(str) + 1;
  loop_setting = new char[n];
  strcpy(loop_setting,str);

  str = (char *) "gaussian";
  n = strlen(str) + 1;
  dist_setting = new char[n];
  strcpy(dist_setting,str);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"min") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal prd command");
      etol = atof(arg[iarg+1]);
      ftol = atof(arg[iarg+2]);
      maxiter = atoi(arg[iarg+3]);
      maxeval = atoi(arg[iarg+4]);
      if (maxiter < 0) error->all(FLERR,"Illegal prd command");
      iarg += 5;

    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal prd command");
      temp_flag = 1;
      temp_dephase = atof(arg[iarg+1]);
      if (temp_dephase <= 0.0) error->all(FLERR,"Illegal prd command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal prd command");
      delete [] loop_setting;
      delete [] dist_setting;

      if (strcmp(arg[iarg+1],"all") == 0) loop_setting = NULL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop_setting = NULL;
      else if (strcmp(arg[iarg+1],"geom") == 0) loop_setting = NULL;
      else error->all(FLERR,"Illegal prd command");
      int n = strlen(arg[iarg+1]) + 1;
      loop_setting = new char[n];
      strcpy(loop_setting,arg[iarg+1]);

      if (strcmp(arg[iarg+2],"uniform") == 0) dist_setting = NULL;
      else if (strcmp(arg[iarg+2],"gaussian") == 0) dist_setting = NULL;
      else error->all(FLERR,"Illegal prd command");
      n = strlen(arg[iarg+2]) + 1;
      dist_setting = new char[n];
      strcpy(dist_setting,arg[iarg+2]);

      iarg += 3;
    } else error->all(FLERR,"Illegal prd command");
  }
}
