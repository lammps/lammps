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
CommandStyle(mdi/engine, MDIEngine);
// clang-format on
#else

#ifndef LMP_MDI_ENGINE_H
#define LMP_MDI_ENGINE_H

#include "command.h"

namespace LAMMPS_NS {

class MDIEngine : public Command {
 public:
  MDIEngine(LAMMPS *lmp) : Command(lmp) {}
  virtual ~MDIEngine() {}
  void command(int, char **);

 private:
  class FixMDIEngine *mdi_fix;

  char *mdi_md();
  char *mdi_optg();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
