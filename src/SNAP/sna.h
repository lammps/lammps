/* ----------------------------------------------------------------------
   Authors: Aidan Thompson and Christian Trott, Sandia National Labs, 2012
   Property of Sandia National Labs: Not for External Distribution
------------------------------------------------------------------------- */

#ifndef LMP_SNA_H
#define LMP_SNA_H

#include <complex>
#include "pointers.h"
#include <ctime>

namespace LAMMPS_NS {
struct SNA_LOOPINDICES {
  int j1, j2, j, ma, mb, ma1, ma2, mb1, mb2;
};

struct SNA_LOOPINDICES_J {
  int j1, j2, j;
};

class SNA : protected Pointers {

public:
  SNA(LAMMPS*, double, int, int, int, double, int);

  SNA(LAMMPS* lmp) : Pointers(lmp) {};
  ~SNA();
  void build_indexlist();
  void init();
  void test();
  double memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  void compute_ui(int);
  void compute_ui_omp(int, int);
  void compute_zi();
  void compute_zi_omp(int);
  void compute_bi();
  void copy_bi2bvec();

  // functions for derivatives

  void compute_duidrj(double*, double, double);
  void compute_dbidrj();
  void compute_dbidrj_nonsymm();
  void copy_dbi2dbvec();
  double compute_sfac(double, double);
  double compute_dsfac(double, double);

  double* timers;
  timespec starttime, endtime;
  int print;
  int counter;

  //per sna class instance for OMP use

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
  double*** uarraytot_r_b, *** uarraytot_i_b;
  double***** zarray_r_b, ***** zarray_i_b;
  double*** uarray_r, *** uarray_i;

private:
  double rmin0, rfac0;

  //use indexlist instead of loops, constructor generates these
  SNA_LOOPINDICES* idx;
  int idx_max;
  SNA_LOOPINDICES_J* idxj;
  int idxj_max;
  // data for bispectrum coefficients

  double***** cgarray;
  double** rootpqarray;
  double*** barray;

  // derivatives of data

  double**** duarray_r, **** duarray_i;
  double**** dbarray;

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
  void add_uarraytot_omp(double, double, double);
  void compute_uarray(double, double, double,
                      double, double);
  void compute_uarray_omp(double, double, double,
                          double, double, int);
  double factorial(int);
  double deltacg(int, int, int);
  int compute_ncoeff();
  void compute_duarray(double, double, double,
                       double, double, double, double, double);

  // if number of atoms are small use per atom arrays 
  // for twojmax arrays, rij, inside, bvec
  // this will increase the memory footprint considerably, 
  // but allows parallel filling and reuse of these arrays
  int use_shared_arrays;

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine
  int switch_flag;

  // Self-weight
  double wself;

};

}

#endif
