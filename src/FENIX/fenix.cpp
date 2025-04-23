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
  input_world = world;

  fenix::args::FenixInitArgs fenix_args;
  fenix_args.in_comm = input_world;
  fenix_args.out_comm = &resilient_world;
  fenix_args.spares = spare_ranks;

  // Only non-spare ranks leave init
  fenix::init(fenix_args);

  world = resilient_world;
  MPI_Comm_size(world, &comm->nprocs);

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
      if (fenix::role() != fenix::INITIAL_RANK) {
        input->one("jump " + restart_file + " " + restart_label);
      }
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
    } else {
      error->all(FLERR, "Invalid argument fenix {}", arg[i]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::fault_handler(){
  world = resilient_world;
  MPI_Comm_rank(world, &comm->me);
  int me = comm->me;

  if(comm->me == 0){
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
      if (me == 0) {
        error->warning(FLERR,
          "not enough spare ranks to maintain initial communicator size"
        );
      }
    } else error->all(FLERR, "Fenix recovery error code: {}", status);
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::recover(){
  while(true) {
    try {
      recover_impl();
      break;
    } catch(const fenix::CommException& e) {}
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::recover_impl(){
  // Destroy and recreate data structures
  lmp->destroy();
  delete input;

  int me;
  MPI_Comm_rank(world, &me);

  if(me == 0 && !infile){
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

}
