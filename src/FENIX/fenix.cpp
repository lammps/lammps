/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "FENIX/fenix.h"
#include "comm.h"
#include "universe.h"
#include "error.h"
#include "utils.h"
#include "input.h"
#include "output.h"
#include "domain.h"

#include <fenix.hpp>

#include <signal.h>
#include <cstring>
#include <cassert>

namespace LAMMPS_NS {

Fenix* Fenix::active_controller = nullptr;

/* ---------------------------------------------------------------------- */

Fenix::Fenix(LAMMPS* lmp) : Command(lmp) {
  spare_ranks = 0;
  resilient_world = MPI_COMM_NULL;
  input_world = MPI_COMM_NULL;
  restart_file = "SELF";
};

/* ---------------------------------------------------------------------- */

void Fenix::command(int narg, char** arg) {
  if (!active_controller) {
    active_controller = new Fenix(lmp);
    active_controller->parse_args(narg, arg);
    active_controller->init();
  } else {
    // If we reach the Fenix command a second time, it's an opportunity
    // to update args but is otherwise a no-op
    active_controller->parse_args(narg, arg);
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::init() {
  assert(active_controller == this);
  if(domain->box_exist){
    error->all(FLERR, "Cannot set up Fenix after simulation box is defined");
  }

  if (!universal) {
    input_world = world;
  } else {
    spare_ranks = universe->procs_per_world[universe->nworlds-1];
    input_world = universe->uworld;
    universe->nworlds--;
  }

  fenix::args::FenixInitArgs fenix_args;
  fenix_args.in_comm = input_world;
  fenix_args.out_comm = &resilient_world;
  fenix_args.spares = spare_ranks;

  // Only non-spare ranks leave init
  fenix::init(fenix_args);

  if (!universal) {
    world = resilient_world;
    MPI_Comm_size(world, &comm->nprocs);
  } else {
    universe->uworld = resilient_world;
    MPI_Comm_size(universe->uworld, &universe->nprocs);
  }

  // Callback is invoked after rebuilding a comm and before an exception
  // is thrown, which lets us log and detect recovery failure before any
  // custom exception handling is invoked
  fenix::callback_register([&](MPI_Comm, int){
      fault_handler();
  });

  // Recovered ranks (spares that just replaced a failed rank) manually
  // invoke the handler and enter recovery
  if(fenix::role() == fenix::RECOVERED_RANK){
    fenix::callback_invoke_all();
    recover();
  }

  // Hijack the input processing to wrap it in an exception handler
  while(true){
    try {
      if (fenix::role() == fenix::INITIAL_RANK) try_setup_universe();
      else input->one("jump " + restart_file + " " + restart_label);
      input->file();
      Fenix_Finalize();
      break;
    } catch(const fenix::CommException& e) {
      recover();
    }
  }

  assert(active_controller == this);
  active_controller = nullptr;
  delete this;
}

/* ---------------------------------------------------------------------- */

void Fenix::parse_args(int narg, char** arg){
  spare_ranks = 0;
  for(int i = 0; i < narg; i++){
    if(strcmp(arg[i], "spares") == 0) {
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix spares", error);
      spare_ranks = utils::inumeric(FLERR, arg[++i], false, lmp);
    } else if(strcmp(arg[i], "restart_file") == 0) {
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix restart_file", error);
      restart_file = std::string(arg[++i]);
    } else if(strcmp(arg[i], "restart_label") == 0) {
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix restart_label", error);
      restart_label = std::string(arg[++i]);
    } else if(strcmp(arg[i], "universal") == 0) {
      if(input_world != MPI_COMM_NULL && !universal) {
        error->all(FLERR, "Cannot make fenix universal after initial setup.");
      }
      universal = true;
    } else {
      error->all(FLERR, "Invalid argument fenix {}", arg[i]);
    }
  }
  if (universal && spare_ranks) error->all(FLERR,
    "When running in universal mode, spares are taken from the final world and "
    "cannot be specified as a number."
  );
}

/* ---------------------------------------------------------------------- */

void Fenix::fault_handler(){
  int me;
  if (!universal) {
    world = resilient_world;
    MPI_Comm_rank(world, &comm->me);
    me = comm->me;
  } else {
    universe->uworld = resilient_world;
    MPI_Comm_rank(universe->uworld, &universe->me);
    me = universe->me;
  }

  if(me == 0){
    std::string fail_list = "[";
    for(auto i : fenix::fail_list()) fail_list += std::to_string(i) + ",";
    fail_list = fail_list.substr(0, fail_list.size()-1) + "]";

    utils::logmesg(lmp,
        "\n\n\nFenix recovering from rank failure(s) {} with {} spare ranks"
        " remaining\n\n", fail_list, fenix::nspare()
    );
  }

  int status = fenix::error();
  if(status != FENIX_SUCCESS){
    if(status == FENIX_WARNING_SPARE_RANKS_DEPLETED){
      // Shrinking supported for non-universal, but not currently for universal.
      auto msg = "not enough spare ranks to maintain initial communicator size";
      if (universal) error->universe_all(FLERR, msg);
      else if(me == 0) error->warning(FLERR, msg);
    } else {
      auto msg = "Fenix recovery error code: " + std::to_string(status);
      if (universal) error->universe_all(FLERR, msg);
      else error->all(FLERR, msg);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::recover(){
  while(true) {
    try {
      try_recover();
      break;
    } catch(const fenix::CommException& e) {}
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::try_recover(){
  // Destroy and recreate data structures
  lmp->destroy();
  delete input;

  if (universal) try_setup_universe();

  int me;
  MPI_Comm_rank(world, &me);

  if(infile) fclose(infile);
  if(screen && universe->existflag) fclose(screen);
  if(logfile && universe->existflag) fclose(logfile);
  if(universe->uscreen) fclose(universe->uscreen);
  if(universe->ulogfile) fclose(universe->ulogfile);
  infile = screen = logfile = universe->uscreen = universe->ulogfile = nullptr;

  if(me == 0){
    // New rank 0 needs to reconfigure I/O args
    std::string lg, scrn, ulg, uscrn;
    for(int i = 0; i < lmp->num_in_arg; i++){
      if(!strcmp("-in",lmp->in_args[i])) infile = fopen(lmp->in_args[++i], "r");
      else if(!strcmp("-plog",lmp->in_args[i]))    lg = lmp->in_args[++i];
      else if(!strcmp("-pscreen",lmp->in_args[i])) scrn = lmp->in_args[++i];
      else if(!strcmp("-log",lmp->in_args[i]))     ulg = lmp->in_args[++i];
      else if(!strcmp("-screen",lmp->in_args[i]))  uscrn = lmp->in_args[++i];
    }
    if(!infile) infile = stdin;
    if(ulg.empty()) ulg = "log.lammps";
    if(lg.empty()) lg = ulg;
    if(scrn.empty()) scrn = uscrn.empty() ? "screen" : uscrn;

    int ume = universe->existflag ? universe->me : me;
    if(ume == 0 && uscrn.empty()) universe->uscreen = stdout;
    else if(ume == 0 && uscrn == "none") universe->uscreen = nullptr;
    else if(ume == 0) universe->uscreen = fopen(uscrn.c_str(), "a");
    if(ume == 0 && ulg == "none") universe->ulogfile = nullptr;
    else universe->ulogfile = fopen(ulg.c_str(), "a");
    if(!universe->existflag){
      screen = universe->uscreen;
      logfile = universe->ulogfile;
    } else {
      int uid = universe->iworld;
      if(scrn == "none") screen = nullptr;
      else screen = fopen(fmt::format("{}.{}", scrn, uid).c_str(), "a");
      if(lg == "none") logfile = nullptr;
      else logfile = fopen(fmt::format("{}.{}", lg, uid).c_str(), "a");
    }
  }

  input = new Input(lmp, lmp->num_in_arg, lmp->in_args);
  lmp->create();
  lmp->post_create();
}

void Fenix::try_setup_universe(){
  // Set mpi_error_uniform to construct, so communicator construction functions
  // always return a uniform error code. This doesn't technically do anything
  // yet, but it will make communicator construction safe even when the
  // construction function is interrupted by a failure once support is added in
  // MPI implementations.
  MPI_Info comm_info;
  MPI_Comm_get_info(universe->uworld, &comm_info);
  MPI_Info_set(comm_info, "mpi_error_uniform", "construct");
  MPI_Info_free(&comm_info);

  // Make sure universe->uni2orig remains accurate
  MPI_Group g_uorig, g_uworld;
  MPI_Comm_group(universe->uorig, &g_uorig);
  MPI_Comm_group(universe->uworld, &g_uworld);
  int* uworld_ranks = new int[universe->nprocs];
  for(int i = 0; i < universe->nprocs; i++) uworld_ranks[i] = i;
  MPI_Group_translate_ranks(
    g_uworld, universe->nprocs, uworld_ranks, g_uorig, universe->uni2orig
  );
  delete[] uworld_ranks;
  MPI_Group_free(&g_uorig);
  MPI_Group_free(&g_uworld);

  // Now rebuild each world's comm
  if (world != MPI_COMM_NULL) MPI_Comm_free(&world);
  int me = universe->me, iworld = -1, count = 0;
  for(int i = 0; i < universe->nworlds; i++){
    count += universe->procs_per_world[i];
    if (me < count) {
      iworld = i;
      break;
    }
  }
  // Only spares should have moved iworlds, this is tested and a proper printout
  // is made in the fault handler.
  assert(iworld == universe->iworld || universe->iworld == universe->nworlds);
  universe->iworld = iworld;
  int err = MPI_Comm_split(universe->uworld, universe->iworld, 0, &world);
  if (err != MPI_SUCCESS) {
    // Needed until "mpi_error_uniform" is supported - after that we can just
    // return (the split will be retried until success)
    error->universe_one(FLERR, "Could not safely recreate universe worlds");
  }
}

}
