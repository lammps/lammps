/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_ORIENT_FCC_H
#define FIX_ORIENT_FCC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOrientFCC : public Fix {
 public:
  struct Nbr {              // neighbor info for each owned and ghost atom
    int n;                  // # of closest neighbors (up to 12)
    int id[12];             // IDs of each neighbor
                            // if center atom is owned, these are local IDs
                            // if center atom is ghost, these are global IDs
    double xismooth[12];    // distance weighting factor for each neighbors
    double dxi[12][3];      // d order-parameter / dx for each neighbor
    double duxi;            // d Energy / d order-parameter for atom
  };

  struct Sort {             // data structure for sorting to find 12 closest 
    int id;                 // ID of neighbor atom
    double rsq;             // distance between center and neighbor atom
    double delta[3];        // displacement between center and neighbor atom
    double xismooth;        // distance weighting factor
  };

  FixOrientFCC(class LAMMPS *, int, char **);
  ~FixOrientFCC();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  double compute_scalar();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  double memory_usage();

 private:
  int me;
  double PI;
  int nlevels_respa;

  int direction_of_motion;         // 1 = center shrinks, 0 = center grows
  int nstats;                      // stats output every this many steps
  double a;                        // lattice parameter
  double Vxi;                      // potential value
  double uxif_low;                 // cut-off fraction, low order parameter
  double uxif_high;                // cut-off fraction, high order parameter
  char *xifilename, *chifilename;  // file names for 2 crystal orientations

  bool use_xismooth;
  double Rxi[12][3],Rchi[12][3],half_xi_chi_vec[2][6][3];
  double xiid,xi0,xi1,xicutoffsq,cutsq,added_energy,total_added_e;
  int half_fcc_nn,nmax;

  Nbr *nbr;
  Sort *sort;
  class NeighList *list;

  void find_best_ref(double *, int, double &, double *);
  static int compare(const void *, const void *);
};

}

#endif
