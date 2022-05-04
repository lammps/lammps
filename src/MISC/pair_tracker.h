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
  void transfer_history(double *, double *) override;

 protected:
  int sizeflag;
  int history;
  int size_history;
  int neighprev;
  double **cut;
  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;
  int freeze_group_bit;
  int store_local_freq;

  char *id_fix_store_local;
  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;
  class FixStoreLocal *fix_store_local;

  int **type_filter;
  double tmin;

  int nvalues, ncount;
  double *output_data;
  typedef void (PairTracker::*FnPtrPack)(int, int, int, double *);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_id1(int, int, int, double *);
  void pack_id2(int, int, int, double *);
  void pack_time_created(int, int, int, double *);
  void pack_time_broken(int, int, int, double *);
  void pack_time_total(int, int, int, double *);
  void pack_x(int, int, int, double *);
  void pack_y(int, int, int, double *);
  void pack_z(int, int, int, double *);
  void pack_rmin(int, int, int, double *);
  void pack_rave(int, int, int, double *);

  void process_data(int, int, double *);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
