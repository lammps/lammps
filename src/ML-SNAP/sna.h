/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
  SNA(LAMMPS *, double, int, double, int, int, int, int, int, int, int);

  SNA(LAMMPS *lmp) : Pointers(lmp){};
  ~SNA() override;
  void build_indexlist();
  void init();
  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  void compute_ui(int, int);
  void compute_zi();
  void compute_yi(const double *);
  void compute_yterm(int, int, int, const double *);
  void compute_bi(int);

  // functions for derivatives

  void compute_duidrj(int);
  void compute_dbidrj();
  void compute_deidrj(double *);
  double compute_sfac(double, double, double, double);
  double compute_dsfac(double, double, double, double);

  // public bispectrum data

  int twojmax;
  double *blist;
  double **dblist;

  // short neighbor list data

  void grow_rij(int);
  int nmax;    // allocated size of short lists

  double **rij;      // short rij list
  int *inside;       // short neighbor list
  double *wj;        // short weight list
  double *rcutij;    // short cutoff list

  // only allocated for switch_inner_flag=1

  double *sinnerij;    // short inner cutoff midpoint list
  double *dinnerij;    // short inner half-width list

  // only allocated for chem_flag=1

  int *element;    // short element list [0,nelements)

 private:
  double rmin0, rfac0;

  // data for bispectrum coefficients

  SNA_ZINDICES *idxz;
  SNA_BINDICES *idxb;

  double **rootpqarray;
  double *cglist;
  int ***idxcg_block;

  double *ulisttot_r, *ulisttot_i;
  double **ulist_r_ij, **ulist_i_ij;    // short u list
  int *idxu_block;

  double *zlist_r, *zlist_i;
  int ***idxz_block;

  int ***idxb_block;

  double **dulist_r, **dulist_i;
  int elem_duarray;    // element of j in derivative

  double *ylist_r, *ylist_i;
  int idxcg_max, idxu_max, idxz_max, idxb_max;

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void print_clebsch_gordan();
  void init_rootpqarray();
  void zero_uarraytot(int);
  void add_uarraytot(double, int);
  void compute_uarray(double, double, double, double, double, int);
  double deltacg(int, int, int);
  void compute_ncoeff();
  void compute_duarray(double, double, double, double, double, double, double, double, int);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Sets the style for the inner switching function
  // 0 = none
  // 1 = cosine
  int switch_inner_flag;

  // Self-weight
  double wself;

  int bzero_flag;       // 1 if bzero subtracted from barray
  double *bzero;        // array of B values for isolated atoms
  int bnorm_flag;       // 1 if barray divided by j+1
  int chem_flag;        // 1 for multi-element bispectrum components
  int wselfall_flag;    // 1 for adding wself to all element labelings
  int nelements;        // number of elements
  int ndoubles;         // number of multi-element pairs
  int ntriples;         // number of multi-element triplets
};

}    // namespace LAMMPS_NS

#endif
