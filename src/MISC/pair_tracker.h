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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tracker,PairTracker);
// clang-format on
#else

#ifndef LMP_PAIR_TRACKER_H
#define LMP_PAIR_TRACKER_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTracker : public Pair {
 public:
  PairTracker(class LAMMPS *);
  ~PairTracker() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  double atom2cut(int) override;
  double radii2cut(double, double) override;

 protected:
  int sizeflag;
  int history;
  int size_history;
  int neighprev;
  double **cut;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;
  int freeze_group_bit;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;
  class FixPairTracker *fix_pair_tracker;

  void transfer_history(double *, double *) override;
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair tracker requires atom attribute radius for finite cutoffs

The atom style defined does not have these attributes.

E: Could not find pair fix neigh history ID

The associated fix neigh/history is missing

E: Cannot use pair tracker without fix pair/tracker

This pairstyle requires one to define a pair/tracker fix

*/
