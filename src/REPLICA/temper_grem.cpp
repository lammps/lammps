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
   Contributing author: David Stelter (BU)
------------------------------------------------------------------------- */

#include "temper_grem.h"

#include "compute.h"
#include "domain.h"
#include "error.h"
#include "finish.h"
#include "fix.h"
#include "fix_grem.h"
#include "fix_nh.h"
#include "force.h"
#include "integrate.h"
#include "modify.h"
#include "random_park.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define TEMPER_DEBUG 0

/* ---------------------------------------------------------------------- */

TemperGrem::TemperGrem(LAMMPS *lmp) : Command(lmp) {}

/* ---------------------------------------------------------------------- */

TemperGrem::~TemperGrem()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete[] set_lambda;
  delete[] lambda2world;
  delete[] world2lambda;
  delete[] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void TemperGrem::command(int narg, char **arg)
{
  if (universe->nworlds == 1)
    error->universe_all(FLERR,"More than one processor partition required for temper/grem command");
  if (domain->box_exist == 0)
    error->universe_all(FLERR,"Temper/grem command before simulation box is defined");
  if (narg != 7 && narg != 8) error->universe_all(FLERR,"Illegal temper/grem command");

  int nsteps = utils::inumeric(FLERR,arg[0],false,lmp);
  nevery = utils::inumeric(FLERR,arg[1],false,lmp);
  double lambda = utils::numeric(FLERR,arg[2],false,lmp);

  // ignore temper command, if walltime limit was already reached

  if (timer->is_timeout()) return;

  // Get and check if gREM fix exists and is correct style

  auto ifix = modify->get_fix_by_id(arg[3]);
  if (!ifix) error->universe_all(FLERR,fmt::format("Tempering fix ID {} is not defined", arg[3]));

  fix_grem = dynamic_cast<FixGrem*>(ifix);
  if (!fix_grem || (strcmp(ifix->style,"grem") != 0))
    error->universe_all(FLERR,"Tempering temperature fix is of incorrect style");

  // Check input values lambdas should be equal, assign other gREM values
  if (lambda != fix_grem->lambda)
    error->universe_all(FLERR,"Lambda from tempering and fix in the same world must be the same");
  double eta = fix_grem->eta;
  double h0 = fix_grem->h0;
  double pressref = 0;

  // Get and check for nh fix

  FixNH *nh = dynamic_cast<FixNH *>(modify->get_fix_by_id(arg[4]));
  if (!nh)
    error->universe_all(FLERR,fmt::format("Fix {} for Nose-Hoover fix does not exist", arg[4]));

  // get result from nvt vs npt check from fix_grem

  int pressflag = fix_grem->pressflag;

  // fix_grem does all the checking...

  if (pressflag) {
    int dummy;
    auto p_start = (double *) nh->extract("p_start",dummy);
    pressref = p_start[0];
  }

  seed_swap = utils::inumeric(FLERR,arg[5],false,lmp);
  seed_boltz = utils::inumeric(FLERR,arg[6],false,lmp);

  my_set_lambda = universe->iworld;
  if (narg == 8) my_set_lambda = utils::inumeric(FLERR,arg[7],false,lmp);
  if ((my_set_lambda < 0) || (my_set_lambda >= universe->nworlds))
    error->universe_one(FLERR,"Illegal temperature index");

  // swap frequency must evenly divide total # of timesteps

  if (nevery <= 0)
    error->universe_all(FLERR,"Invalid frequency in temper command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in temper command");

  // setup for long tempering run

  update->whichflag = 1;
  timer->init_timeout();

  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0) error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage

  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap

  Compute *pe_compute = modify->get_compute_by_id("thermo_pe");
  if (!pe_compute) error->all(FLERR,"Tempering could not find thermo_pe compute");

  pe_compute->addstep(update->ntimestep + nevery);

  // create MPI communicator for root proc from each world

  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG

  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = nullptr;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i

  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set lambdas
  // allgather tempering arg "lambda" across root procs
  // bcast from each root to other procs in world

  set_lambda = new double[nworlds];
  if (me == 0) MPI_Allgather(&lambda,1,MPI_DOUBLE,set_lambda,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_lambda,nworlds,MPI_DOUBLE,0,world);

  // create world2lambda only on root procs from my_set_lambda
  // create lambda2world on root procs from world2lambda,
  //   then bcast to all procs within world

  world2lambda = new int[nworlds];
  lambda2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_lambda,1,MPI_INT,world2lambda,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) lambda2world[world2lambda[i]] = i;
  }
  MPI_Bcast(lambda2world,nworlds,MPI_INT,0,world);

  // if restarting tempering, reset lambda target of Fix to current my_set_lambda

  if (narg == 8) {
    double new_lambda = set_lambda[my_set_lambda];
    fix_grem->lambda = new_lambda;
  }

  // setup tempering runs

  int i,which,partner,swap,partner_set_lambda,partner_world;
  double pe,weight,weight_partner,weight_cross, weight_cross_partner;
  double volume,enth,new_lambda,boltz_factor;

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up tempering ...\n");

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }

  timer->init();
  timer->barrier_start();

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps

    timer->init_timeout();
    update->integrate->run(nevery);

    // check for timeout across all procs

    int my_timeout=0;
    int any_timeout=0;
    if (timer->is_timeout()) my_timeout=1;
    MPI_Allreduce(&my_timeout, &any_timeout, 1, MPI_INT, MPI_SUM, universe->uworld);
    if (any_timeout) {
      timer->force_timeout();
      break;
    }

    // compute PE
    // notify compute it will be called at next swap

    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + nevery);

    // which = which of 2 kinds of swaps to do (0,1)

    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_lambda = which set lambda I am partnering with for this swap

    if (which == 0) {
      if (my_set_lambda % 2 == 0) partner_set_lambda = my_set_lambda + 1;
      else partner_set_lambda = my_set_lambda - 1;
    } else {
      if (my_set_lambda % 2 == 1) partner_set_lambda = my_set_lambda + 1;
      else partner_set_lambda = my_set_lambda - 1;
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps

    partner = -1;
    if (me == 0 && partner_set_lambda >= 0 && partner_set_lambda < nworlds) {
      partner_world = lambda2world[partner_set_lambda];
      partner = world2root[partner_world];
    }

    // swap with a partner, only root procs in each world participate
    // hi proc sends PE to low proc
    // lo proc make Boltzmann decision on whether to swap
    // lo proc communicates decision back to hi proc

    swap = 0;
    if (partner != -1) {
      // compute weights
      volume = domain->xprd * domain->yprd * domain->zprd;
      enth = pe + (pressref * volume);
      weight = log(set_lambda[my_set_lambda] + (eta*(enth - h0)));
      weight_cross = log(set_lambda[partner_set_lambda] + (eta*(enth - h0)));

      if (me_universe > partner) {
        MPI_Send(&weight,1,MPI_DOUBLE,partner,0,universe->uworld);
        MPI_Send(&weight_cross,1,MPI_DOUBLE,partner,0,universe->uworld);
      }
      else {
        MPI_Recv(&weight_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
        MPI_Recv(&weight_cross_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
      }

      if (me_universe < partner) {
        boltz_factor = (weight + weight_partner - weight_cross - weight_cross_partner) *
            (1 / (boltz * eta));
        if (boltz_factor >= 0.0) swap = 1;
        else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
      }

      if (me_universe < partner)
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);

#if TEMPER_DEBUG
      if (me_universe < partner)
        fprintf(universe->uscreen,"SWAP %d & %d: yes = %d,Ts = %d %d, PEs = %g %g, Bz = %g %g\n",
                me_universe,partner,swap,my_set_lambda,partner_set_lambda,
                weight,weight_partner,boltz_factor,exp(boltz_factor));
#endif
    }

    // bcast swap result to other procs in my world

    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // if my world swapped, all procs in world reset temp target of Fix

    if (swap) {
      new_lambda = set_lambda[partner_set_lambda];
      fix_grem->lambda = new_lambda;
    }

    // update my_set_lambda and lambda2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world

    if (swap) my_set_lambda = partner_set_lambda;
    if (me == 0) {
      MPI_Allgather(&my_set_lambda,1,MPI_INT,world2lambda,1,MPI_INT,roots);
      for (i = 0; i < nworlds; i++) lambda2world[world2lambda[i]] = i;
    }
    MPI_Bcast(lambda2world,nworlds,MPI_INT,0,world);

    // print out current swap status

    if (me_universe == 0) print_status();
  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void TemperGrem::print_status()
{
  std::string status = std::to_string(update->ntimestep);
  for (int i = 0; i < nworlds; i++)
    status += " " + std::to_string(world2lambda[i]);

  status += "\n";

  if (universe->uscreen) fputs(status.c_str(), universe->uscreen);
  if (universe->ulogfile) {
    fputs(status.c_str(), universe->ulogfile);
    fflush(universe->ulogfile);
  }
}
