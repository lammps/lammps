/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

  // static method for MDI to callback to
  // when LAMMPS is a driver which launches a plugin engine

  static int plugin_wrapper(void *, MDI_Comm, void *);
};

}    // namespace LAMMPS_NS

#endif
