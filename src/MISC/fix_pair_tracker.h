/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(pair/tracker,FixPairTracker);
// clang-format on
#else

#ifndef LMP_FIX_PAIR_TRACKING_H
#define LMP_FIX_PAIR_TRACKING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPairTracker : public Fix {
 public:
  FixPairTracker(class LAMMPS *, int, char **);
  ~FixPairTracker();
  int setmask();
  void init();
  void post_force(int);
  double memory_usage();
  void lost_contact(int, int, double, double, double, double);

 private:
  int nvalues, nmax;
  int index_i, index_j;
  double tmin, rmin, rsum, time_initial, nstep_initial;

  double *vector;
  double **array;
  int **type_filter;

  int ncount;

  void reallocate(int);

  typedef void (FixPairTracker::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_id1(int);
  void pack_id2(int);

  void pack_time_created(int);
  void pack_time_broken(int);
  void pack_time_total(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);

  void pack_rmin(int);
  void pack_rave(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid keyword in fix pair/tracker command

Self-explanatory.

*/
