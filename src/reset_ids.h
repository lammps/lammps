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

CommandStyle(reset_ids,ResetIDs)

#else

#ifndef LMP_RESET_IDS_H
#define LMP_RESET_IDS_H

#include "pointers.h"

namespace LAMMPS_NS {

class ResetIDs : protected Pointers {
 public:
  ResetIDs(class LAMMPS *);
  void command(int, char **);

 private:
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Reset_ids command before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

UNDOCUMENTED

E: Cannot use reset_ids unless atoms have IDs

UNDOCUMENTED

E: Reset_ids missing %d bond topology atom IDs - use comm_modify cutoff

UNDOCUMENTED

*/
