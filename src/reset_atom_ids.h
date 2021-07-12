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
CommandStyle(reset_atom_ids,ResetIDs);
// clang-format on
#else

#ifndef LMP_RESET_IDS_H
#define LMP_RESET_IDS_H

#include "command.h"

namespace LAMMPS_NS {

class ResetIDs : public Command {
 public:
  struct AtomRvous {
    bigint ibin;
    int proc, ilocal;
    double x[3];
  };

  struct IDRvous {
    tagint newID;
    int ilocal;
  };

#if defined(LMP_QSORT)
  // static variable across all ResetID objects, for qsort callback
  static AtomRvous *sortrvous;
#endif

  ResetIDs(class LAMMPS *);
  void command(int, char **);

 private:
  bigint binlo, binhi;

  // callback functions for rendezvous communication

  static int sort_bins(int, char *, int &, int *&, char *&, void *);

  void sort();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Reset_ids command before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

UNDOCUMENTED

E: Cannot use reset_atom_ids unless atoms have IDs

UNDOCUMENTED

E: Reset_ids missing %d bond topology atom IDs - use comm_modify cutoff

UNDOCUMENTED

*/
