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

#include "fenix.hpp"

#include <signal.h>
#include <fmt/ranges.h>

namespace LAMMPS_NS {

Fenix* Fenix::active_controller = nullptr;

/* ---------------------------------------------------------------------- */

Fenix::Fenix(LAMMPS* lmp) : Command(lmp), chkpt(lmp) {
  restart_pos = 0;
  chkpt_interval = 0;
  spare_ranks = 0;

  resilient_world = MPI_COMM_NULL;

  orig_input = input;
};

/* ---------------------------------------------------------------------- */

void Fenix::command(int narg, char** arg) {
  if(active_controller){
    if(active_controller->chkpt_interval){
      chkpt.recover();
    } else if(!active_controller->restart_file.empty()){
      input->one("read_restart " + active_controller->restart_file);
    }
    if(!active_controller->jump_cmd.empty()){
      input->one(active_controller->jump_cmd);
    }
    return;
  } else {
    active_controller = this;
  }

  parse_args(narg, arg);

  MPI_Comm full_world = world;

  fenix::args::FenixInitArgs fenix_args;
  fenix_args.in_comm = full_world;
  fenix_args.out_comm = &resilient_world;
  fenix_args.spares = spare_ranks;

  // Only non-spare ranks leave init
  fenix::init(fenix_args);

  world = resilient_world;
  MPI_Comm_size(world, &comm->nprocs);

  // Automatically repair LAMMPS on faults.
  fenix::callback_register([&](MPI_Comm, int){
      fault_handler();
  });

  if(fenix::role() == fenix::INITIAL_RANK){
    setup_checkpointing();
  } else if(fenix::role() == fenix::RECOVERED_RANK){
    fenix::callback_invoke_all();
  }

  // Hijack the input processing to wrap it in an exception handler
  while(true){
    try{
      input->file();
      Fenix_Finalize();
      break;
    } catch(const fenix::CommException& e) {
      // Error handling is done within the callback registered to fenix
    }
  }

  if(output->restart == &chkpt) output->restart = nullptr;
  if(input != orig_input){
    delete input;
    input = orig_input;
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::parse_args(int narg, char** arg){
  for(int i = 0; i < narg; i++){
    if(strcmp(arg[i], "spares") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix spares", error);
      spare_ranks = utils::inumeric(FLERR, arg[++i], false, lmp);
    } else if(strcmp(arg[i], "checkpoint_every") == 0) {
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix checkpoint_every", error);
      chkpt_interval = utils::inumeric(FLERR, arg[++i], false, lmp);
    } else if(strcmp(arg[i], "restart_here") == 0) {
      //Warning: don't use this yet! Probably more feasible for localized recovery
      if(infile == stdin) error->one(FLERR,
        "fenix restart_here invalid when infile is stdin"
      );
      if(comm->me == 0) restart_pos = platform::ftell(infile);
      else restart_pos = -1;
    } else if(strcmp(arg[i], "restart_jump") == 0) {
      if(i+2 >= narg) utils::missing_cmd_args(FLERR, "fenix restart_jump", error);
      jump_cmd  = "jump " + std::string(arg[++i]);
      jump_cmd += " " + std::string(arg[++i]);
    } else if(strcmp(arg[i], "restart_file") == 0) {
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fenix restart_file", error);
      restart_file = arg[++i];
    } else {
      error->all(FLERR, "Invalid argument fenix {}", arg[i]);
    }
  }

  if(chkpt_interval != 0 && !restart_file.empty()) error->all(FLERR,
    "Invalid arguments fenix checkpoint_every and restart_file are incompatible"
  );
}

/* ---------------------------------------------------------------------- */

void Fenix::fault_handler(){
  world = resilient_world;
  MPI_Comm_rank(world, &comm->me);
  int me = comm->me;

  if(comm->me == 0){
    utils::logmesg(lmp,
        "\n\n\nFenix recovering from rank failure(s) [{}] with {} spare ranks"
        " remaining\n\n",
        fmt::join(fenix::fail_list(), ", "), fenix::nspare()
    );
  }

  int status = fenix::error();
  if(status != FENIX_SUCCESS){
    if(status == FENIX_WARNING_SPARE_RANKS_DEPLETED){
      error->all(FLERR,
        "Fenix recovery failed, not enough spare ranks to maintain initial"
        " communicator size"
      );
    } else error->all(FLERR, "Fenix recovery error code: {}", status);
  }

  // Remove chkpt from output so it doesn't delete it when recreated
  if(output->restart == &chkpt)
    output->restart = nullptr;

  // Destroy and recreate data structures
  lmp->destroy();

  if(input != orig_input){
    delete input;
  }

  if(me == 0 && !infile){
    // New rank 0 needs to reconfigure I/O args
    std::string lg = "", scrn = "", ulg = "", uscrn = "";
    for(int i = 0; i < lmp->num_in_arg; i++){
      if(!strcmp("-in",lmp->in_args[i])) infile = fopen(lmp->in_args[++i], "r");
      else if(!strcmp("-plog",lmp->in_args[i]))    lg = lmp->in_args[++i];
      else if(!strcmp("-pscreen",lmp->in_args[i])) scrn = lmp->in_args[++i];
      else if(!strcmp("-log",lmp->in_args[i]))     ulg = lmp->in_args[++i];
      else if(!strcmp("-screen",lmp->in_args[i]))  uscrn = lmp->in_args[++i];
    }
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

  // Must be after recreating lmp
  setup_checkpointing();

  if(chkpt_interval && restart_pos != 0){
    chkpt.restore_metadata();
  }

  if(infile && infile != stdin){
    if(restart_pos == -1) error->one(FLERR,
        "Fenix was unable to recover the input file's restart position "
        "after a failure on rank 0"
    );

    platform::fseek(infile, restart_pos);
  }

  // If restarting from beginning, we'll restore LAMMPS when we get back to
  // the fenix command. Otherwise we restore now and continue from after the
  // fenix command.
  if(restart_pos != 0){
    chkpt.recover();
  }
}

/* ---------------------------------------------------------------------- */

void Fenix::setup_checkpointing(){
  if(!chkpt_interval) return;
  chkpt.create_group();

  // We'll just piggyback on the existing checkpoint framework
  // Note that this breaks all kinds of stuff if the user later uses the
  // builtin restart command, since Output will delete chkpt
  output->restart_flag = 1;
  output->restart_flag_single = 1;
  output->restart_every_single = chkpt_interval;
  delete[] output->restart1;
  output->restart1 = new char[1];
  output->restart1[0] = '\0';

  if(output->restart != &chkpt){
    delete output->restart;
    output->restart = &chkpt;
  }
}

}
