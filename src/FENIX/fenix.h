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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(fenix,Fenix);
// clang-format on
#else

#ifndef LMP_FENIX_H
#define LMP_FENIX_H

#include "command.h"
#include "FENIX/fenix_checkpoint.h"

namespace LAMMPS_NS {

class Fenix : public Command {
 public:
  Fenix(class LAMMPS *);
  void command(int, char**) override;

  void parse_args(int, char**);

  void fault_handler();

  void setup_checkpointing();

  int spare_ranks;
  int chkpt_interval;
  std::string jump_cmd;

  MPI_Comm resilient_world;

  Input* orig_input;

  //Metadata is checkpointed
  bigint metadata[1];
  bigint& restart_pos = metadata[0];

  FenixCheckpoint chkpt;

  static Fenix* active_controller;
};

}

#endif
#endif
