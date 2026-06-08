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

namespace LAMMPS_NS {

class Fenix : public Command {
 public:
  Fenix(class LAMMPS *);
  void command(int, char**) override;

  void parse_args(int, char**);
  void init();

  void fault_handler();

  void recover();
  void try_recover();

  void try_setup_universe();

  int spare_ranks;
  bool universal;
  std::string restart_file;
  std::string restart_label;

  MPI_Comm input_world;
  MPI_Comm resilient_world;

  static Fenix* active_controller;
};

}

#endif
#endif
