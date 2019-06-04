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

#include <complex>
#include "pointers.h"
#include <ctime>

namespace LAMMPS_NS {

struct SNA_ZINDICES {
  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju;
  double betaj;
};

struct SNA_BINDICES {
  int j1, j2, j;
};

class SNA : protected Pointers {

public:
  SNA(LAMMPS*, double, int, int, double, int, int);

  SNA(LAMMPS* lmp) : Pointers(lmp) {};
  ~SNA();
  void build_indexlist();
  void init();
  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  void compute_ui(int);
  void compute_zi();
  void compute_yi(const double*);
  void compute_bi();
  void copy_bi2bvec();

  // functions for derivatives

  void compute_duidrj(double*, double, double);
  void compute_dbidrj();
  void compute_deidrj(double*);
  void copy_dbi2dbvec();
  double compute_sfac(double, double);
  double compute_dsfac(double, double);

  double* bvec, ** dbvec;
  double** rij;
  int* inside;
  double* wj;
  double* rcutij;
  int nmax;

  void grow_rij(int);

  int twojmax, diagonalstyle;
  double*** uarraytot_r, *** uarraytot_i;
  double***** zarray_r, ***** zarray_i;
  double*** yarray_r, *** yarray_i;
  double*** uarray_r, *** uarray_i;

private:
  double rmin0, rfac0;

  // use indexlist instead of loops, constructor generates these

  SNA_ZINDICES* idxz;
  SNA_BINDICES* idxb;
  int idxcg_max, idxu_max, idxz_max, idxb_max;

  // data for bispectrum coefficients

  double***** cgarray;
  double** rootpqarray;
  double*** barray;

  // derivatives of data

  double**** duarray_r, **** duarray_i;
  double**** dbarray;

  double* cglist;  
  int*** idxcg_block; 

  double* ulisttot_r, * ulisttot_i;
  double* ulist_r, * ulist_i;
  int* idxu_block;

  double* zlist_r, * zlist_i;
  int*** idxz_block;

  int*** idxb_block;

  double** dulist_r, ** dulist_i;
  double* ylist_r, * ylist_i;

  static const int nmaxfactorial = 167;
  static const double nfac_table[];
  double factorial(int);

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void init_rootpqarray();
  void jtostr(char*, int);
  void mtostr(char*, int, int);
  void print_clebsch_gordan(FILE*);
  void zero_uarraytot();
  void addself_uarraytot(double);
  void add_uarraytot(double, double, double);
  void compute_uarray(double, double, double,
                      double, double);
  double deltacg(int, int, int);
  int compute_ncoeff();
  void compute_duarray(double, double, double,
                       double, double, double, double, double);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Self-weight
  double wself;

  int bzero_flag; // 1 if bzero subtracted from barray
  double *bzero;  // array of B values for isolated atoms
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid argument to factorial %d

N must be >= 0 and <= 167, otherwise the factorial result is too
large.

*/
