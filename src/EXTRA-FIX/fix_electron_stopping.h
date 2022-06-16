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
/* ----------------------------------------------------------------------
   Electronic stopping power
   Contributing authors: K. Avchaciov and T. Metspalu
   Information: k.avchachov@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(electron/stopping,FixElectronStopping);
// clang-format on
#else

#ifndef LMP_FIX_ELECTRON_STOPPING_H
#define LMP_FIX_ELECTRON_STOPPING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixElectronStopping : public Fix {
 public:
  FixElectronStopping(class LAMMPS *, int, char **);
  ~FixElectronStopping() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;
  void init_list(int, class NeighList *) override;
  double compute_scalar() override;

 private:
  void read_table(const char *);
  void grow_table();

  double Ecut;                  // cutoff energy
  double SeLoss, SeLoss_all;    // electronic energy loss
  int SeLoss_sync_flag;         // sync done since last change?

  int maxlines;              // max number of lines in table
  int table_entries;         // number of table entries actually read
  double **elstop_ranges;    // [ 0][i]: energies
                             // [>0][i]: stopping powers per type

  char *idregion;          // region id
  class Region *region;    // region pointer if used, else NULL
  int minneigh;            // minimum number of neighbors

  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
