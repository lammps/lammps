/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(nn/snap,PairNNSNAP)

#else

#ifndef LMP_PAIR_NN_SNAP_H
#define LMP_PAIR_NN_SNAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNNSNAP : public Pair {
public:
  PairNNSNAP(class LAMMPS *);
  ~PairNNSNAP();
  virtual void compute(int, int);
  void compute_regular(int, int);
  void compute_optimized(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();

  double rcutfac, quadraticflag; // declared public to workaround gcc 4.9
  int ncoeff;                    //  compiler bug, manifest in KOKKOS package

protected:
  int ncoeffq, ncoeffall;
  double **bvec, ***dbvec;
  class SNA** sna;
  int nmax;
  int nthreads;
  virtual void allocate();
  void read_files(char *, char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);
  double extra_cutoff();
  void load_balance();
  void set_sna_to_shared(int snaid,int i);
  void build_per_atom_arrays();

  void compute_beta();

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

#ifdef TIMING_INFO
  //  timespec starttime, endtime;
  double timers[4];
#endif

  double rcutmax;               // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  double *radelem;              // element radii
  double *wjelem;               // elements weights
  double **coeffelem;           // element bispectrum coefficients
  double** beta;                // betas for all atoms in list
  int *map;                     // mapping from atom types to elements
  int twojmax, diagonalstyle, switchflag, bzeroflag;
  double rfac0, rmin0, wj1, wj2;
  int rcutfacflag, twojmaxflag; // flags for required parameters
  int beta_max;                 // length of beta
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Communication cutoff too small for SNAP micro load balancing

This can happen if you change the neighbor skin after your pair_style
command or if your box dimensions grow during a run. You can set the
cutoff explicitly via the comm_modify cutoff command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Must set number of threads via package omp command

Because you are using the USER-OMP package, set the number of threads
via its settings, not by the pair_style snap nthreads setting.

W: Communication cutoff is too small for SNAP micro load balancing, increased to %lf

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Incorrect SNAP coeff file

UNDOCUMENTED

E: Incorrect SNAP parameter file

The file cannot be parsed correctly, check its internal syntax.

E: Pair style SNAP requires newton pair on

See the newton command.  This is a restriction to use the SNAP
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open SNAP coefficient file %s

The specified SNAP coefficient file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in SNAP coefficient file

Incorrect number of words per line in the coefficient file.

E: Cannot open SNAP parameter file %s

The specified SNAP parameter file cannot be opened.  Check that the
path and name are correct.

E: Incorrect format in SNAP parameter file

Incorrect number of words per line in the parameter file.

E: Did not find all elements in SNAP coefficient file.

One or more elements listed in the pair_coeff command were not found in the coefficient file.

*/
