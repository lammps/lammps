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
CommandStyle(mdi/plugin,MDIPlugin);
// clang-format on
#else

#ifndef LMP_MDI_PLUGIN_H
#define LMP_MDI_PLUGIN_H

#include "command.h"
#include "mdi.h"

namespace LAMMPS_NS {

class MDIPlugin : public Command {
 public:
  MDIPlugin(LAMMPS *lmp) : Command(lmp) {}
  void command(int, char **) override;

 private:
  char *lammps_command;
  class Fix *fixptr;

  static int plugin_wrapper(void *, MDI_Comm, void *);
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
