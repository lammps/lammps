// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "tad.h"

#include "atom.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "fix_event_tad.h"
#include "fix_store_atom.h"
#include "force.h"
#include "integrate.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "neb.h"
#include "neighbor.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

TAD::TAD(LAMMPS *lmp) : Command(lmp)
{
  deltconf = deltstop = deltfirst = 0.0;
}

/* ---------------------------------------------------------------------- */

TAD::~TAD()
{
  memory->sfree(fix_event_list);
  if (neb_logfilename != nullptr) delete[] neb_logfilename;
  delete[] min_style;
  delete[] min_style_neb;
}

/* ----------------------------------------------------------------------
   perform TAD simulation on root proc
   other procs only used for NEB calcs
------------------------------------------------------------------------- */

void TAD::command(int narg, char **arg)
{
  fix_event_list = nullptr;
  n_event_list = 0;
  nmax_event_list = 0;
  nmin_event_list = 10;

  // error checks

  if (domain->box_exist == 0)
    error->all(FLERR,"Tad command before simulation box is defined");
  if (universe->nworlds == 1)
    error->all(FLERR,"Cannot use TAD with a single replica for NEB");
  if (universe->nworlds != universe->nprocs)
    error->all(FLERR,"Can only use TAD with 1-processor replicas for NEB");
  if (atom->sortfreq > 0)
    error->all(FLERR,"Cannot use TAD with atom_modify sort enabled for NEB");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Cannot use TAD unless atom map exists for NEB");

  if (narg < 7) error->universe_all(FLERR,"Illegal tad command");

  nsteps = utils::inumeric(FLERR,arg[0],false,lmp);
  t_event = utils::inumeric(FLERR,arg[1],false,lmp);
  templo = utils::numeric(FLERR,arg[2],false,lmp);
  temphi = utils::numeric(FLERR,arg[3],false,lmp);
  delta_conf = utils::numeric(FLERR,arg[4],false,lmp);
  tmax = utils::numeric(FLERR,arg[5],false,lmp);

  char *id_compute = utils::strdup(arg[6]);

  // quench minimizer is set by min_style command
  // NEB minimizer is set by options, default = quickmin

  min_style = utils::strdup(update->minimize_style);

  options(narg-7,&arg[7]);

  // total # of timesteps must be multiple of t_event

  if (t_event <= 0) error->universe_all(FLERR,"Invalid t_event in tad command");
  if (nsteps % t_event)
    error->universe_all(FLERR,"TAD nsteps must be multiple of t_event");

  if (delta_conf <= 0.0 || delta_conf >= 1.0)
    error->universe_all(FLERR,"Invalid delta_conf in tad command");

  if (tmax <= 0.0)
    error->universe_all(FLERR,"Invalid tmax in tad command");

  // deltconf = (ln(1/delta))/freq_min (timestep units)

  deltconf = -log(delta_conf)*tmax/update->dt;

  // local storage

  int me_universe = universe->me;
  int nprocs_universe = universe->nprocs;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  delta_beta = (1.0/templo - 1.0/temphi) / force->boltz;
  ratio_beta = templo/temphi;

  // create FixEventTAD object to store last event

  fix_event = dynamic_cast<FixEventTAD *>(modify->add_fix("tad_event all EVENT/TAD"));

  // create FixStoreAtom object to store revert state

  fix_revert = dynamic_cast<FixStoreAtom *>(modify->add_fix("tad_revert all STORE/ATOM 7 0 0 0"));

  // create Finish for timing output

  finish = new Finish(lmp);

  // assign FixEventTAD to event-detection compute
  // necessary so it will know atom coords at last event

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID for TAD");
  compute_event = modify->compute[icompute];
  compute_event->reset_extra_compute_fix("tad_event");

  // reset reneighboring criteria since will perform minimizations

  neigh_every = neighbor->every;
  neigh_delay = neighbor->delay;
  neigh_dist_check = neighbor->dist_check;

  if (neigh_every != 1 || neigh_delay != 0 || neigh_dist_check != 1) {
    if (me_universe == 0)
      error->warning(FLERR,"Resetting reneighboring criteria during TAD");
  }

  neighbor->every = 1;
  neighbor->delay = 0;
  neighbor->dist_check = 1;

  // initialize TAD as if one long dynamics run

  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  update->restrict_output = 1;
  if (update->laststep < 0)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // set minimize style for quench

  char *args[1];
  args[0] = min_style;

  update->create_minimize(1,args,1);

  // init minimizer settings and minimizer itself

  update->etol = etol;
  update->ftol = ftol;
  update->max_eval = maxeval;

  update->minimize->init();

  // perform TAD simulation

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up TAD ...\n");

  if (me_universe == 0) {
    if (universe->uscreen)
      fprintf(universe->uscreen,
              "Step CPU N M Status Barrier Margin t_lo delt_lo\n");
    if (universe->ulogfile)
      fprintf(universe->ulogfile,
              "Step CPU N M Status Barrier Margin t_lo delt_lo\n");
  }

  ulogfile_lammps = universe->ulogfile;
  uscreen_lammps = universe->uscreen;
  ulogfile_neb = nullptr;
  uscreen_neb = nullptr;
  if (me_universe == 0 && neb_logfilename)
    ulogfile_neb = fopen(neb_logfilename,"w");

  // store hot state and quenched event, only on replica 0

  // need this line if quench() does only setup_minimal()
  // update->minimize->setup();

  // this should work with if statement uncommented, but does not
  // if (universe->iworld == 0) {

  fix_event->store_state_quench();
  quench();

  timer->init();
  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);
  fix_event->store_event_tad(update->ntimestep);
  log_event(0);
  fix_event->restore_state_quench();

  // do full init/setup

  update->whichflag = 1;
  lmp->init();
  update->integrate->setup(1);

  // main loop: look for events until out of time
  // (1) dynamics, store state, quench, check event, restore state
  // (2) if event, perform NEB, record in fix_event_list
  // (3) if confident, pick earliest event

  nbuild = ndanger = 0;
  time_neb = time_dynamics = time_quench = time_comm = time_output = 0.0;

  timer->init();
  timer->barrier_start();
  time_start = timer->get_wall(Timer::TOTAL);

  int confident_flag, event_flag;

  if (universe->iworld == 0) {
    while (update->ntimestep < update->endstep) {

      // initialize list of possible events

      initialize_event_list();
      confident_flag = 0;

      while (update->ntimestep < update->endstep) {
        event_flag = 0;
        while (update->ntimestep < update->endstep) {

          dynamics();
          fix_event->store_state_quench();
          quench();

          event_flag = check_event();
          MPI_Bcast(&event_flag,1,MPI_INT,0,universe->uworld);

          if (event_flag) break;

          // restore hot state

          fix_event->restore_state_quench();

          // store hot state in revert

          store_state();
        }
        if (!event_flag) break;

        add_event();

        perform_neb(n_event_list-1);
        compute_tlo(n_event_list-1);
        confident_flag = check_confidence();
        MPI_Bcast(&confident_flag,1,MPI_INT,0,universe->uworld);
        if (confident_flag) break;
        if (universe->iworld == 0) revert_state();
      }
      if (!confident_flag) break;

      perform_event(event_first);

      // need to sync timestep with TAD

      MPI_Bcast(&(update->ntimestep),1,MPI_INT,0,universe->uworld);

      int restart_flag = 0;
      if (output->restart_flag && universe->iworld == 0) {
        if (output->restart_every_single &&
            fix_event->event_number % output->restart_every_single == 0)
          restart_flag = 1;
        if (output->restart_every_double &&
            fix_event->event_number % output->restart_every_double == 0)
          restart_flag = 1;
      }

      // full init/setup since are starting after event

      update->whichflag = 1;
      lmp->init();
      update->integrate->setup(1);

    // write restart file of hot coords

      if (restart_flag) {
        timer->barrier_start();
        output->write_restart(update->ntimestep);
        timer->barrier_stop();
        time_output += timer->get_wall(Timer::TOTAL);
      }
    }

  } else {

    while (update->ntimestep < update->endstep) {
      confident_flag = 0;
      while (update->ntimestep < update->endstep) {
        event_flag = 0;
        while (update->ntimestep < update->endstep) {
          update->ntimestep += t_event;
          MPI_Bcast(&event_flag,1,MPI_INT,0,universe->uworld);

          if (event_flag) break;
        }
        if (!event_flag) break;
        perform_neb(-1);
        MPI_Bcast(&confident_flag,1,MPI_INT,0,universe->uworld);
        if (confident_flag) break;
      }
      if (!confident_flag) break;

      // need to sync timestep with TAD

      MPI_Bcast(&(update->ntimestep),1,MPI_INT,0,universe->uworld);
    }
  }

  // set total timers and counters so Finish() will process them

  timer->set_wall(Timer::TOTAL, time_start);
  timer->barrier_stop();

  timer->set_wall(Timer::NEB, time_neb);
  timer->set_wall(Timer::DYNAMICS, time_dynamics);
  timer->set_wall(Timer::QUENCH, time_quench);
  timer->set_wall(Timer::REPCOMM, time_comm);
  timer->set_wall(Timer::REPOUT, time_output);

  neighbor->ncalls = nbuild;
  neighbor->ndanger = ndanger;

  if (me_universe == 0) {
    auto mesg = fmt::format("Loop time of {} on {} procs for {} steps with {} atoms\n",
                            timer->get_wall(Timer::TOTAL), nprocs_universe, nsteps,atom->natoms);
    if (universe->uscreen) fmt::print(universe->uscreen, mesg);
    if (universe->ulogfile) fmt::print(universe->ulogfile, mesg);
  }

  if ((me_universe == 0) && ulogfile_neb) fclose(ulogfile_neb);

  if (me == 0) utils::logmesg(lmp,"\nTAD done\n");

  finish->end(3);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
  update->restrict_output = 0;

  // reset reneighboring criteria

  neighbor->every = neigh_every;
  neighbor->delay = neigh_delay;
  neighbor->dist_check = neigh_dist_check;


  delete[] id_compute;
  delete finish;
  modify->delete_fix("tad_event");
  modify->delete_fix("tad_revert");
  delete_event_list();

  compute_event->reset_extra_compute_fix(nullptr);
}

/* ----------------------------------------------------------------------
   single short dynamics run
------------------------------------------------------------------------- */

void TAD::dynamics()
{
  update->whichflag = 1;
  update->nsteps = t_event;

  lmp->init();
  update->integrate->setup(1);
  // this may be needed if don't do full init
  //modify->addstep_compute_all(update->ntimestep);
  int ncalls = neighbor->ncalls;

  timer->barrier_start();
  update->integrate->run(t_event);
  timer->barrier_stop();
  time_dynamics += timer->get_wall(Timer::TOTAL);

  nbuild += neighbor->ncalls - ncalls;
  ndanger += neighbor->ndanger;

  update->integrate->cleanup();
  finish->end(0);
}

/* ----------------------------------------------------------------------
   quench minimization
------------------------------------------------------------------------- */

void TAD::quench()
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

  // full init works

  lmp->init();
  update->minimize->setup();

  // partial init does not work

  //modify->addstep_compute_all(update->ntimestep);
  //update->minimize->setup_minimal(1);

  int ncalls = neighbor->ncalls;

  timer->barrier_start();
  update->minimize->run(maxiter);
  timer->barrier_stop();
  time_quench += timer->get_wall(Timer::TOTAL);

  if (neighbor->ncalls == ncalls) quench_reneighbor = 0;
  else quench_reneighbor = 1;

  update->minimize->cleanup();
  finish->end(1);

  // reset timestep as if quench did not occur
  // clear timestep storage from computes, since now invalid

  update->ntimestep = ntimestep_hold;
  update->endstep = update->laststep = endstep_hold;
  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->clearstep();
}

/* ----------------------------------------------------------------------
   check for an event
   return 0 if no event
   return 1 if event
------------------------------------------------------------------------- */

int TAD::check_event()
{
  int flag;

  flag = 0;
  if (compute_event->compute_scalar() > 0.0) flag = 1;

  return flag;
}

/* ----------------------------------------------------------------------
   universe proc 0 prints event info
------------------------------------------------------------------------- */

void TAD::log_event(int ievent)
{
  timer->set_wall(Timer::TOTAL, time_start);
  if (universe->me == 0) {
    double tfrac = 0.0;
    auto mesg = fmt::format("{} {:.3f} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f}\n",
                            fix_event->event_timestep, timer->elapsed(Timer::TOTAL),
                            fix_event->event_number, ievent, "E ", fix_event->ebarrier,
                            tfrac, fix_event->tlo, deltfirst);

    if (universe->uscreen) fmt::print(universe->uscreen, mesg);
    if (universe->ulogfile) fmt::print(universe->ulogfile, mesg);
  }

  // dump snapshot of quenched coords
  // must reneighbor and compute forces before dumping
  // addstep_compute_all ensures eng/virial are calculated if needed

  if (output->ndump && universe->iworld == 0) {
    timer->barrier_start();
    modify->addstep_compute_all(update->ntimestep);
    update->integrate->setup_minimal(1);
    // must reset whichflag so that computes won't fail.
    update->whichflag = 1;
    output->write_dump(update->ntimestep);
    update->whichflag = 0;
    timer->barrier_stop();
    time_output += timer->get_wall(Timer::TOTAL);
  }

}

/* ----------------------------------------------------------------------
   parse optional parameters at end of TAD input line
------------------------------------------------------------------------- */

void TAD::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal tad command");

  // set defaults

  etol = 0.1;
  ftol = 0.1;
  maxiter = 40;
  maxeval = 50;

  etol_neb = 0.01;
  ftol_neb = 0.01;
  n1steps_neb = 100;
  n2steps_neb = 100;
  nevery_neb = 10;

  min_style_neb = utils::strdup("quickmin");
  dt_neb = update->dt;
  neb_logfilename = nullptr;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"min") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal tad command");
      etol = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ftol = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      maxiter = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      maxeval = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      if (maxiter < 0 || maxeval < 0 ||
          etol < 0.0 || ftol < 0.0 )
        error->all(FLERR,"Illegal tad command");
      iarg += 5;

    } else if (strcmp(arg[iarg],"neb") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal tad command");
      etol_neb = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ftol_neb = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      n1steps_neb = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      n2steps_neb = utils::inumeric(FLERR,arg[iarg+4],false,lmp);
      nevery_neb = utils::inumeric(FLERR,arg[iarg+5],false,lmp);
      if (etol_neb < 0.0 || ftol_neb < 0.0 ||
          n1steps_neb < 0 || n2steps_neb < 0 ||
          nevery_neb < 0) error->all(FLERR,"Illegal tad command");
      iarg += 6;

    } else if (strcmp(arg[iarg],"neb_style") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal tad command");
      delete[] min_style_neb;
      min_style_neb = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"neb_step") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal tad command");
      dt_neb = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (dt_neb <= 0.0) error->all(FLERR,"Illegal tad command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"neb_log") == 0) {
      delete[] neb_logfilename;
      if (iarg+2 > narg) error->all(FLERR,"Illegal tad command");
      if (strcmp(arg[iarg+1],"none") == 0) neb_logfilename = nullptr;
      else {
        neb_logfilename = utils::strdup(arg[iarg+1]);
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal tad command");
  }
}

/* ----------------------------------------------------------------------
   perform NEB calculation
------------------------------------------------------------------------- */

void TAD::perform_neb(int ievent)
{

  double **x = atom->x;
  int nlocal = atom->nlocal;

  double *buf_final;
  memory->create(buf_final,3*nlocal,"tad:buffinal");

  // set system to quenched state of event ievent

  if (universe->iworld == 0) {

    fix_event_list[ievent]->restore_event();

    int ii = 0;
    for (int i = 0; i < nlocal; i++) {
      buf_final[ii++] = x[i][0];
      buf_final[ii++] = x[i][1];
      buf_final[ii++] = x[i][2];
    }
  }

  MPI_Bcast(buf_final,3*nlocal,MPI_DOUBLE,universe->root_proc[0],
            universe->uworld);

  double *buf_init;
  memory->create(buf_init,3*nlocal,"tad:bufinit");

  // set system to quenched state of fix_event

  if (universe->iworld == 0) {

    fix_event->restore_event();

    int ii = 0;
    for (int i = 0; i < nlocal; i++) {
      buf_init[ii++] = x[i][0];
      buf_init[ii++] = x[i][1];
      buf_init[ii++] = x[i][2];
    }
  }

  MPI_Bcast(buf_init,3*nlocal,MPI_DOUBLE,universe->root_proc[0],
            universe->uworld);

  // create FixNEB object to support NEB

  fix_neb = (Fix *) modify->add_fix("neb all neb 1.0");

  // switch minimize style to quickmin for NEB

  char *args[1];
  args[0] = min_style_neb;
  update->create_minimize(1,args,1);

  // create NEB object

  neb = new NEB(lmp,etol_neb,ftol_neb,n1steps_neb,
                n2steps_neb,nevery_neb,buf_init,buf_final);

  // free up temporary arrays

  memory->destroy(buf_init);
  memory->destroy(buf_final);

  // run NEB with NEB timestep

  int beginstep_hold = update->beginstep;
  int endstep_hold = update->endstep;
  int ntimestep_hold = update->ntimestep;
  int nsteps_hold = update->nsteps;

  if (universe->me == 0) {
    universe->ulogfile = ulogfile_neb;
    universe->uscreen = uscreen_neb;
  }

  // had to bypass timer interface
  // because timer->array is reset inside neb->run()

  //    timer->barrier_start();
  //    neb->run();
  //    timer->barrier_stop();
  //    time_neb += timer->get_wall(Timer::TOTAL);

  MPI_Barrier(world);
  double time_tmp = platform::walltime();

  double dt_hold = update->dt;
  update->dt = dt_neb;
  neb->run();
  update->dt = dt_hold;

  MPI_Barrier(world);
  time_neb += platform::walltime() - time_tmp;

  if (universe->me == 0) {
    universe->ulogfile = ulogfile_lammps;
    universe->uscreen = uscreen_lammps;
  }

  // extract barrier energy from NEB

  if (universe->iworld == 0)
    fix_event_list[ievent]->ebarrier = neb->ebf;

  update->beginstep = update->firststep = beginstep_hold;
  update->endstep = update->laststep = endstep_hold;
  update->ntimestep = ntimestep_hold;
  update->nsteps = nsteps_hold;

  // switch minimize style back for quench

  args[0] = min_style;
  update->create_minimize(1,args,1);

  update->etol = etol;
  update->ftol = ftol;

  // clean up

  modify->delete_fix("neb");
  delete neb;
}

/* ----------------------------------------------------------------------
   check if confidence criterion for tstop is satisfied
   return 0 if not satisfied
   return 1 if satisfied
------------------------------------------------------------------------- */

int TAD::check_confidence()
{
  int flag;

  // update stopping time

  deltstop = deltconf*pow(deltfirst/deltconf, ratio_beta);

  flag = 0;
  if (deltstop < update->ntimestep - fix_event->event_timestep) flag = 1;

  return flag;
}

/* ----------------------------------------------------------------------
   store state in fix_revert
------------------------------------------------------------------------- */

void TAD::store_state()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double **astore = fix_revert->astore;

  for (int i = 0; i < nlocal; i++) {
    astore[i][0] = x[i][0];
    astore[i][1] = x[i][1];
    astore[i][2] = x[i][2];
    astore[i][3] = v[i][0];
    astore[i][4] = v[i][1];
    astore[i][5] = v[i][2];
    *((imageint *) &astore[i][6]) = image[i];
  }
}

/* ----------------------------------------------------------------------
   restore state archived in fix_revert
   flip sign of velocities to reflect back to starting state
------------------------------------------------------------------------- */

void TAD::revert_state()
{
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  double **astore = fix_revert->astore;

  for (int i = 0; i < nlocal; i++) {
    x[i][0] = astore[i][0];
    x[i][1] = astore[i][1];
    x[i][2] = astore[i][2];
    v[i][0] = -astore[i][3];
    v[i][1] = -astore[i][4];
    v[i][2] = -astore[i][5];
    image[i] = *((imageint *) &astore[i][6]);
  }
}

/* ----------------------------------------------------------------------
   initialize list of possible events
------------------------------------------------------------------------- */

void TAD::initialize_event_list() {

  // First delete old events, if any

  delete_event_list();

  // Create new list

  n_event_list = 0;
  grow_event_list(nmin_event_list);
}

/* ----------------------------------------------------------------------
   delete list of possible events
------------------------------------------------------------------------- */

void TAD::delete_event_list() {

  for (int i = 0; i < n_event_list; i++)
    modify->delete_fix(fmt::format("tad_event_{}",i));

  memory->sfree(fix_event_list);
  fix_event_list = nullptr;
  n_event_list = 0;
  nmax_event_list = 0;

}

/* ----------------------------------------------------------------------
   add event
------------------------------------------------------------------------- */

void TAD::add_event()
{
  if (n_event_list == nmax_event_list)
    grow_event_list(nmax_event_list+nmin_event_list);

  // create FixEventTAD object to store possible event

  int ievent = n_event_list++;
  fix_event_list[ievent]
    = dynamic_cast<FixEventTAD *>(modify->add_fix(fmt::format("tad_event_{} all EVENT/TAD", ievent)));

  // store quenched state for new event

  fix_event_list[ievent]->store_event_tad(update->ntimestep);

  // store hot state for new event

  fix_event->restore_state_quench();
  fix_event_list[ievent]->store_state_quench();
}

/* ----------------------------------------------------------------------
   compute cold time for event ievent
------------------------------------------------------------------------- */

void TAD::compute_tlo(int ievent)
{
  double deltlo,delthi,ebarrier;

  ebarrier = fix_event_list[ievent]->ebarrier;
  delthi = fix_event_list[ievent]->event_timestep
    - fix_event->event_timestep;
  deltlo = delthi*exp(ebarrier*delta_beta);
  fix_event_list[ievent]->tlo = fix_event->tlo + deltlo;

  // update first event

  auto  statstr = (char *) "D ";

  if (ievent == 0) {
    deltfirst = deltlo;
    event_first = ievent;
    statstr = (char *) "DF";
  } else if (deltlo < deltfirst) {
    deltfirst = deltlo;
    event_first = ievent;
    statstr = (char *) "DF";
  }

  // first-replica output about each event

  timer->set_wall(Timer::TOTAL, time_start);
  if (universe->me == 0) {
    double tfrac = 0.0;
    if (ievent > 0) tfrac = delthi/deltstop;
    auto mesg = fmt::format("{} {:.3f} {} {} {} {:.3f} {:.3f} {:.3f} {:.3f}\n",
                            fix_event_list[ievent]->event_timestep, timer->elapsed(Timer::TOTAL),
                            fix_event->event_number, ievent, statstr, ebarrier, tfrac,
                            fix_event->tlo, deltlo);

    if (universe->uscreen) fmt::print(universe->uscreen, mesg);
    if (universe->ulogfile) fmt::print(universe->ulogfile, mesg);
  }
}

/* ----------------------------------------------------------------------
   perform event
------------------------------------------------------------------------- */

void TAD::perform_event(int ievent)
{
  // reset timestep to that of event

  update->ntimestep = fix_event_list[ievent]->event_timestep;

  // Copy event to current event
  // Should really use copy constructor for this
  fix_event->tlo = fix_event_list[ievent]->tlo;
  fix_event->ebarrier = fix_event_list[ievent]->ebarrier;
  fix_event->event_number++;
  fix_event->event_timestep = update->ntimestep;
  fix_event_list[ievent]->restore_event();
  fix_event->store_event_tad(fix_event_list[ievent]->event_timestep);

  // output stats and dump for quenched state

  log_event(ievent);

  // load and store hot state

  fix_event_list[ievent]->restore_state_quench();
  fix_event->store_state_quench();
}

/* ----------------------------------------------------------------------
   Allocate list of pointers to events
------------------------------------------------------------------------- */

void TAD::grow_event_list(int nmax) {
  if (nmax_event_list > nmax) return;
  fix_event_list = (FixEventTAD **)
    memory->srealloc(fix_event_list,nmax*sizeof(FixEventTAD *),"tad:eventlist");
  nmax_event_list = nmax;
}
