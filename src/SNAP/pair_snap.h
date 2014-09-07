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

#ifdef PAIR_CLASS

PairStyle(snap,PairSNAP)

#else

#ifndef LMP_PAIR_SNAP_H
#define LMP_PAIR_SNAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSNAP : public Pair {
public:
  PairSNAP(class LAMMPS *);
  ~PairSNAP();
  void compute(int, int);
  void compute_regular(int, int);
  void compute_optimized(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

protected:
  int ncoeff;
  double **bvec, ***dbvec;
  class SNA** sna;
  int nmax;
  int nthreads;
  void allocate();
  void read_files(char *, char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);
  double extra_cutoff();
  void load_balance();
  void set_sna_to_shared(int snaid,int i);
  void build_per_atom_arrays();

  int schedule_user;
  double schedule_time_guided;
  double schedule_time_dynamic;

  int ncalls_neigh;
  int do_load_balance;
  int ilistmask_max;
  int* ilistmask;
  int ghostinum;
  int ghostilist_max;
  int* ghostilist;
  int ghostnumneigh_max;
  int* ghostnumneigh;
  int* ghostneighs;
  int* ghostfirstneigh;
  int ghostneighs_total;
  int ghostneighs_max;

  int use_optimized;
  int use_shared_arrays;

  int i_max;
  int i_neighmax;
  int i_numpairs;
  int **i_pairs;
  double ***i_rij;
  int **i_inside;
  double **i_wj;
  double **i_rcutij;
  int *i_ninside;
  double ****i_uarraytot_r, ****i_uarraytot_i;
  double ******i_zarray_r, ******i_zarray_i;
  //  timespec starttime, endtime;
  double timers[4];
  double gamma;

  double rcutmax;               // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double *radelem;              // element radii
  double *wjelem;               // elements weights
  double **coeffelem;           // element bispectrum coefficients
  int *map;                     // mapping from atom types to elements
  int twojmax, diagonalstyle, switchflag;
  double rcutfac, rfac0, rmin0, wj1, wj2;
  int rcutfacflag, twojmaxflag; // flags for required parameters
  int gammaoneflag;              // 1 if parameter gamma is 1
};

}

#endif
#endif
