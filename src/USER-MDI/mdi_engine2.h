/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(mdi/engine2, MDIEngine2);
// clang-format on
#else

#ifndef LMP_MDI_ENGINE2_H
#define LMP_MDI_ENGINE2_H

#include "command.h"
#include "mdi.h"

namespace LAMMPS_NS {

class MDIEngine2 : public Command {
 public:
  MDIEngine2(LAMMPS *lmp) : Command(lmp) {}
  virtual ~MDIEngine2() {}
  void command(int, char **);

  int execute_command(const char *command, MDI_Comm mdicomm);
  void engine_node(const char *node);

 private:
  int lmpunits;    // REAL or METAL
  int root;        // 1 for procs 0, otherwise 0
  int mode;        // which mode engine is in (DEFAULT,MD,OPTG,etc)
  char *cmd;
  char *node_driver;   // which node driver is at
  char *node_engine;   // which node engine is at

  MDI_Comm mdicomm;
  class FixMDIEngine2 *mdi_fix;

  bool exit_flag;
  bool local_exit_flag;

  // command to be executed at the target node

  char *target_command;

  char *id_pe;
  char *id_ke;
  class Irregular *irregular;
  class Minimize *minimizer;
  class Compute *pe;
  class Compute *ke;
  double *add_force;

  void mdi_commands();
  void mdi_md();
  void mdi_optg();

  void send_types(Error *);
  void send_labels(Error *);
  void send_masses(Error *);
  void receive_coordinates(Error *);
  void send_coordinates(Error *);
  void send_charges(Error *);
  void send_energy(Error *);
  void send_forces(Error *);
  void send_pe(Error *);
  void send_ke(Error *);
  void receive_forces(Error *, int);
  void send_cell(Error *);
  void receive_cell(Error *);
  void send_celldispl(Error *);
  void receive_celldispl(Error *);
  void exchange_forces();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
