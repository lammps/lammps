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
   Contributing authors: W. Michael Brown, Intel
------------------------------------------------------------------------- */

#ifndef LMP_SNA_INTEL_H
#define LMP_SNA_INTEL_H

#if defined(__AVX512F__)
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)

#include "pointers.h"
#include "intel_buffers.h"
#include "intel_simd.h"

#define SVW 8

#if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
#define SV_for _Pragma("omp simd") _Pragma("vector aligned") for
#else
#define SV_for _Pragma("simd assert") _Pragma("vector aligned") for
#endif
#else
#define SV_for for
#endif

namespace LAMMPS_NS {

struct SNA_ZINDICES {
  int j1, j2, j, ma1min, ma2max, mb1min;
  int mb2max, na, nb, jju;
};

struct SNA_BINDICES {
  int j1, j2, j;
};

#define SNA_DVEC ip_simd::SIMD_double
#define SNA_IVEC ip_simd::SIMD256_int

class SNAIntel : protected Pointers {

 public:
  SNAIntel(LAMMPS *, double, int, double, int, int, int, int, int, int, int);

  SNAIntel(LAMMPS *lmp) : Pointers(lmp){};
  ~SNAIntel() override;
  void build_indexlist();
  void init();
  double memory_usage();

  int ncoeff;

  inline int vector_width() const { return SVW; }

  // functions for bispectrum coefficients

  void compute_ui(const SNA_IVEC &, const SNA_IVEC &, const int max_jnum);
  template <int> void compute_zi_or_yi(const SNA_DVEC *);
  void compute_yi_from_zi(const SNA_DVEC *);
  void compute_yterm(int, int, int, const double *);
  void compute_bi(const SNA_IVEC &);

  // functions for derivatives

  void compute_duidrj(const int, const SNA_IVEC &);
  void compute_deidrj_e(const int, const SNA_IVEC &, SNA_DVEC *);
  void compute_deidrj(const int, const SNA_IVEC &, SNA_DVEC *);
  double compute_sfac(double, double, double, double);
  SNA_DVEC compute_sfac(const SNA_DVEC &, const SNA_DVEC &, const SNA_DVEC &,
                        const SNA_DVEC &);
  inline SNA_DVEC compute_sfac_dsfac(const SNA_DVEC &, const SNA_DVEC &,
                                     const SNA_DVEC &, const SNA_DVEC &,
                                     SNA_DVEC &);

  // public bispectrum data

  int twojmax;
  SNA_DVEC *blist;
  double **dblist;

  // short neighbor list data

  void grow_rij(int);
  int nmax;    // allocated size of short lists

  SNA_DVEC **rij;      // short rij list
  SNA_IVEC *inside;       // short neighbor list
  SNA_DVEC *wj;        // short weight list
  SNA_DVEC *rcutij;    // short cutoff list

  // only allocated for switch_inner_flag=1

  SNA_DVEC *sinnerij;    // short inner cutoff midpoint list
  SNA_DVEC *dinnerij;    // short inner half-width list

  // only allocated for chem_flag=1

  SNA_IVEC *element;    // short element list [0,nelements)

 private:
  double rmin0, rfac0;

  // data for bispectrum coefficients

  SNA_ZINDICES *idxz;
  SNA_BINDICES *idxb;

  double **rootpqarray;
  double *cglist;
  int ***idxcg_block;

  SNA_DVEC *ulisttot_r, *ulisttot_i;
  SNA_DVEC **ulist_r_ij, **ulist_i_ij;    // short u list
  int *idxu_block;

  SNA_DVEC *zlist_r, *zlist_i;
  int ***idxz_block;

  int ***idxb_block;

  SNA_DVEC **dulist_r, **dulist_i;

  SNA_DVEC *ylist_r, *ylist_i;
  int idxcg_max, idxu_max, idxz_max, idxb_max;

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void print_clebsch_gordan();
  void init_rootpqarray();
  void zero_uarraytot(const SNA_IVEC &);
  void add_uarraytot(const SNA_DVEC &, const int, const SNA_IVEC &);
  void compute_uarray(const SNA_DVEC &, const SNA_DVEC &, const SNA_DVEC &,
                      const SNA_DVEC &, const SNA_DVEC &, const int,
                      const SNA_IVEC &);
  double deltacg(int, int, int);
  void compute_ncoeff();
  void compute_duarray(const SNA_DVEC &, const SNA_DVEC &, const SNA_DVEC &,
                       const SNA_DVEC &, const SNA_DVEC &, const SNA_DVEC &,
                       const SNA_DVEC &, const SNA_DVEC &, int,
                       const SNA_IVEC &);
  inline double choose_beta(const int, const int, const int,
                            const int, const int, const int,  int &);

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
#endif

#endif
