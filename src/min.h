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

#ifndef MIN_H
#define MIN_H

#include "pointers.h"

namespace LAMMPS_NS {

class Min : protected Pointers {
 public:
  double einitial,efinal,eprevious;
  double fnorm2_init,fnorminf_init,fnorm2_final,fnorminf_final;
  double alpha_final;
  int niter,neval;
  char *stopstr;

  Min(class LAMMPS *);
  virtual ~Min();
  void init();
  void setup();
  void setup_minimal(int);
  void run(int);
  void cleanup();
  void request(class Pair *, int, double);
  double memory_usage() {return 0.0;}
  void modify_params(int, char **);

  virtual void init_style() {}
  virtual void setup_style() = 0;
  virtual void reset_vectors() = 0;
  virtual int iterate(int) = 0;

 protected:
  int eflag,vflag;            // flags for energy/virial computation
  int virial_style;           // compute virial explicitly or implicitly

  double dmax;                // max dist to move any atom in one step
  int linestyle;              // 0 = backtrack, 1 = quadratic

  int nelist_atom;                  // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_atom;       // list of PE,virial Computes
  class Compute **vlist_global;
  class Compute **vlist_atom;

  int pairflag,torqueflag;
  int triclinic;              // 0 if domain is orthog, 1 if triclinic

  int narray;                       // # of arrays stored by fix_minimize
  class FixMinimize *fix_minimize;  // fix that stores auxiliary data

  class Compute *pe_compute;        // compute for potential energy
  double ecurrent;                  // current potential energy

  double ndoftotal;           // total dof for entire problem

  int n3;                     // local atomic dof
  double *x;                  // variables for atomic dof, as 1d vector
  double *f;                  // force vector for atomic dof, as 1d vector

  int nextra_global;          // # of extra global dof due to fixes
  double *fextra;             // force vector for extra global dof
                              // xextra is stored by fix

  int nextra_atom;            // # of sets of extra per-atom dof
  double **xextra_atom;       // variables for extra per-atom dof sets
  double **fextra_atom;       // force vectors for extra per-atom dof sets
  int *extra_peratom;         // # of per-atom values in each set
  int *extra_nlen;            // total local length of each set, e.g 3*nlocal
  double *extra_max;          // max allowed change in one iter for each set
  class Pair **requestor;     // Pair that requested each set

  int neigh_every,neigh_delay,neigh_dist_check;  // neighboring params

  double energy_force(int);
  void force_clear();

  double compute_force_norm_sqr();
  double compute_force_norm_inf();

  void ev_setup();
  void ev_set(int);

  double fnorm_sqr();
  double fnorm_inf();
};

}

#endif
