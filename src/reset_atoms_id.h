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
CommandStyle(RESET_ATOMS_ID,ResetAtomsID);
// clang-format on
#else

#ifndef LMP_RESET_ATOMS_ID_H
#define LMP_RESET_ATOMS_ID_H

#include "command.h"

namespace LAMMPS_NS {

class ResetAtomsID : public Command {
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

  ResetAtomsID(class LAMMPS *);
  void command(int, char **) override;

 private:
  bigint binlo, binhi;

  // callback functions for rendezvous communication

  static int sort_bins(int, char *, int &, int *&, char *&, void *);

  void sort();
};
}    // namespace LAMMPS_NS
#endif
#endif
