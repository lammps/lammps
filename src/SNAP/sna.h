/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#ifndef LMP_SNA_H
#define LMP_SNA_H

#include "pointers.h"

namespace LAMMPS_NS {

struct SNA_ZINDICES {
  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju;
};

struct SNA_BINDICES {
  int j1, j2, j;
};

class SNA : protected Pointers {

public:
  SNA(LAMMPS*, double, int, double, int, int, int, int, int, int);

  SNA(LAMMPS* lmp) : Pointers(lmp) {};
  ~SNA();
  void build_indexlist();
  void init();
  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  void compute_ui(int, int);
  void compute_zi();
  void compute_yi(const double*);
  void compute_yterm(int, int, int, const double*);
  void compute_bi(int);

  // functions for derivatives

  void compute_duidrj(double*, double, double, int, int);
  void compute_dbidrj();
  void compute_deidrj(double*);
  double compute_sfac(double, double);
  double compute_dsfac(double, double);

  double* blist;
  double** dblist;
  double** rij;
  int* inside;
  double* wj;
  double* rcutij;
  int* element;  // index on [0,nelements)
  int nmax;

  void grow_rij(int);

  int twojmax;

private:
  double rmin0, rfac0;

  // data for bispectrum coefficients

  SNA_ZINDICES* idxz;
  SNA_BINDICES* idxb;
  int idxcg_max, idxu_max, idxz_max, idxb_max;

  double** rootpqarray;
  double* cglist;
  int*** idxcg_block;

  double* ulisttot_r, * ulisttot_i;
  double** ulist_r_ij, ** ulist_i_ij;
  int* idxu_block;

  double* zlist_r, * zlist_i;
  int*** idxz_block;

  int*** idxb_block;

  double** dulist_r, ** dulist_i;
  int elem_duarray; // element of j in derivative
  double* ylist_r, * ylist_i;

  static const int nmaxfactorial = 167;
  static const double nfac_table[];
  double factorial(int);

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void print_clebsch_gordan();
  void init_rootpqarray();
  void zero_uarraytot(int);
  void addself_uarraytot(double, int);
  void add_uarraytot(double, double, double, int, int);
  void compute_uarray(double, double, double,
                      double, double, int);
  double deltacg(int, int, int);
  void compute_ncoeff();
  void compute_duarray(double, double, double,
                       double, double, double, double, double, int);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Self-weight
  double wself;

  int bzero_flag; // 1 if bzero subtracted from barray
  double* bzero;  // array of B values for isolated atoms
  int bnorm_flag; // 1 if barray divided by j+1
  int alloy_flag; // 1 for multi-element bispectrum components
  int wselfall_flag; // 1 for adding wself to all element labelings
  int nelements;  // number of elements
  int ndoubles;   // number of multi-element pairs
  int ntriples;   // number of multi-element triplets
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid argument to factorial %d

N must be >= 0 and <= 167, otherwise the factorial result is too
large.

*/
