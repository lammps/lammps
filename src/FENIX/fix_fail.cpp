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

#include "fix_fail.h"
#include "kill.h"
#include "utils.h"
#include "update.h"
#include "input.h"
#include "comm.h"
#include "error.h"
#include "variable.h"

#include "fenix.hpp"

#include <signal.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixFail::FixFail(LAMMPS* lmp, int narg, char** arg) :
  Fix(lmp, narg, arg)
{
  // Skip 3 args, since they come in as:
  //  arg = <ID> <group-ID> <style> [args...]
  for(int i = 3; i < narg; i++){
    if(strcmp(arg[i], "rank") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fix fail rank", error);
      if(strcmp(arg[++i], "none") == 0){
        fail_rank = -2;
        continue;
      }
      if(strcmp(arg[i], "all") == 0){
        fail_rank = comm->me;
        continue;
      }
      std::string csv = arg[i];

      while(!csv.empty()){
        size_t comma = csv.rfind(",");
        if(comma == std::string::npos) comma = -1;

        size_t dash = csv.find("-", comma+1);
        if(dash == std::string::npos){
          int rank = utils::inumeric(FLERR, &csv[comma+1], false, lmp);
          if(rank == comm->me) fail_rank = rank;
          else if(fail_rank == -1) fail_rank = rank;
        } else {
          int end = utils::inumeric(FLERR, &csv[dash+1], false, lmp);
          csv.erase(dash);
          int begin = utils::inumeric(FLERR, &csv[comma+1], false, lmp);
          if(comm->me >= begin && comm->me <= end) fail_rank = comm->me;
          else if(fail_rank == -1) fail_rank = begin;
        }

        if(comma == -1) csv.clear();
        else csv.erase(comma);
      }
    } else if(strcmp(arg[i], "timestep") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fix fail timestep", error);
      if(strcmp(arg[++i], "none") == 0) fail_timestep = -2;
      else fail_timestep = utils::inumeric(FLERR, arg[i], false, lmp);
    } else if(strcmp(arg[i], "step") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fix fail timestep", error);
      fail_step = utils::fixmask(FLERR, arg[++i], false, lmp);
    } else if(strcmp(arg[i], "var") == 0){
      if(i+2 >= narg) utils::missing_cmd_args(FLERR, "fix fail var", error);
      fail_var = arg[++i];
      fail_var_val = utils::numeric(FLERR, arg[++i], false, lmp);
    } else if(strcmp(arg[i], "wait_only") == 0){
      wait_only = true;
    } else if(strcmp(arg[i], "node") == 0){
      if(i+1 >= narg) utils::missing_cmd_args(FLERR, "fix fail node", error);
      int fail_node = utils::inumeric(FLERR, arg[++i], false, lmp);

      MPI_Comm local_world;
      MPI_Comm_split_type(
        world, MPI_COMM_TYPE_SHARED, comm->me, MPI_INFO_NULL, &local_world
      );
      int local_rank;
      MPI_Comm_rank(local_world, &local_rank);
      MPI_Comm_free(&local_world);

      int is_node_start = local_rank == 0 ? 1 : 0;
      int my_node_index = 0;
      MPI_Scan(&is_node_start, &my_node_index, 1, MPI_INT, MPI_SUM, world);
      my_node_index--;

      if(my_node_index == fail_node) fail_rank = comm->me;
      else fail_rank = -2;
    } else {
      error->all(FLERR, "Invalid argument fix fail {}", arg[i]);
    }
  }

  if(fail_var.empty() && (fail_rank == -1 || fail_timestep == -1)){
    error->all(FLERR, "Invalid arguments, cannot determine where/when to fail");
  }
}

/* ----------------------------------------------------------------------
   Kill this rank if all specified restrictions pass
------------------------------------------------------------------------- */

void FixFail::check_fail(){
  if(skip_double_failure && Fenix::role() == FENIX_ROLE_RECOVERED_RANK) return;
  if(fail_rank != -1 && comm->me != fail_rank) return;
  if(fail_timestep != -1 && update->ntimestep != fail_timestep) return;
  if(!fail_var.empty()){
    double var = input->variable->find(fail_var.c_str());
    if(var < 0) return;
    if(!input->variable->equalstyle(var))
      error->one(FLERR, fmt::format("Variable \"{}\" for fix fail var is invalid style", fail_var));
    if(fail_var_val != input->variable->compute_equal(var)) return;
  }

  if(wait_only){
    //Spin on an MPI barrier until a failure occurs
    while(true) MPI_Barrier(world);
  } else {
    Kill::kill();
  }
}
