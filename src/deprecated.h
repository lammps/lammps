/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(DEPRECATED,Deprecated);
CommandStyle(reset_ids,Deprecated);
CommandStyle(kim_init,Deprecated);
CommandStyle(kim_interactions,Deprecated);
CommandStyle(kim_param,Deprecated);
CommandStyle(kim_property,Deprecated);
CommandStyle(kim_query,Deprecated);
// clang-format on
#else

#ifndef LMP_DEPRECATED_H
#define LMP_DEPRECATED_H

#include "command.h"

namespace LAMMPS_NS {

class Deprecated : public Command {
 public:
  Deprecated(class LAMMPS *lmp) : Command(lmp){};
  void command(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

W: Ignoring unknown or incorrect info command flag

Self-explanatory.  An unknown argument was given to the info command.
Compare your input with the documentation.

E: Unknown name for info package category

Self-explanatory.

E: Unknown name for info newton category

Self-explanatory.

E: Unknown name for info pair category

Self-explanatory.

E: Unknown category for info is_active()

Self-explanatory.

E: Unknown category for info is_available()

Self-explanatory.

E: Unknown category for info is_defined()

Self-explanatory.

*/
