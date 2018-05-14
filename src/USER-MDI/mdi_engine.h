/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(mdi_engine,CommandMDIEngine)

#else

#ifndef LMP_COMMAND_MDI_ENGINE_H
#define LMP_COMMAND_MDI_ENGINE_H

#include "command.h"

namespace LAMMPS_NS {

  class CommandMDIEngine : public Command {
 public:
  CommandMDIEngine(class LAMMPS *);
  virtual ~CommandMDIEngine();
  void command(int, char **);
  int mdi_md();
  int mdi_optg();

private:
  class FixMDIEngine *mdi_fix;
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
