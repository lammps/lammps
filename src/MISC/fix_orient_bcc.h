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

#ifdef FIX_CLASS
// clang-format off
FixStyle(orient/bcc,FixOrientBCC);
// clang-format on
#else

#ifndef LMP_FIX_ORIENT_BCC_H
#define LMP_FIX_ORIENT_BCC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOrientBCC : public Fix {
 public:
  struct Nbr {             // neighbor info for each owned and ghost atom
    int n;                 // # of closest neighbors (up to 8)
    tagint id[8];          // IDs of each neighbor
                           // if center atom is owned, these are local IDs
                           // if center atom is ghost, these are global IDs
    double xismooth[8];    // distance weighting factor for each neighbors
    double dxi[8][3];      // d order-parameter / dx for each neighbor
    double duxi;           // d Energy / d order-parameter for atom
  };

  struct Sort {         // data structure for sorting to find 8 closest
    int id;             // local ID of neighbor atom
    double rsq;         // distance between center and neighbor atom
    double delta[3];    // displacement between center and neighbor atom
    double xismooth;    // distance weighting factor
  };

  FixOrientBCC(class LAMMPS *, int, char **);
  ~FixOrientBCC();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_scalar();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 private:
  int me;
  int ilevel_respa;

  int direction_of_motion;           // 1 = center shrinks, 0 = center grows
  int nstats;                        // stats output every this many steps
  double a;                          // lattice parameter
  double Vxi;                        // potential value
  double uxif_low;                   // cut-off fraction, low order parameter
  double uxif_high;                  // cut-off fraction, high order parameter
  char *xifilename, *chifilename;    // file names for 2 crystal orientations

  bool use_xismooth;
  double Rxi[8][3], Rchi[8][3], half_xi_chi_vec[2][4][3];
  double xiid, xi0, xi1, xicutoffsq, cutsq, added_energy;
  int half_bcc_nn;

  int nmax;          // expose 2 per-atom quantities
  double **order;    // order param and normalized order param

  Nbr *nbr;
  Sort *sort;
  class NeighList *list;

  void find_best_ref(double *, int, double &, double *);
  static int compare(const void *, const void *);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix orient/bcc file open failed

The fix orient/bcc command could not open a specified file.

E: Fix orient/bcc file read failed

The fix orient/bcc command could not read the needed parameters from a
specified file.

E: Fix orient/bcc found self twice

The neighbor lists used by fix orient/bcc are messed up.  If this
error occurs, it is likely a bug, so send an email to the
"developers"_https://www.lammps.org/authors.html.

*/
