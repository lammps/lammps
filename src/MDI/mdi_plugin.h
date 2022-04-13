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

#ifndef LMP_MDI_PLUGIN_H
#define LMP_MDI_PLUGIN_H

#include "mdi.h"
#include "pointers.h"

namespace LAMMPS_NS {

class MDIPlugin : protected Pointers {
 public:
  MDIPlugin(class LAMMPS *, int, char **);

 private:
  char *lammps_command;

  static int plugin_wrapper(void *, MDI_Comm, void *);
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
