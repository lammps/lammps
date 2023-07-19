/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MIN_H
#define LMP_MIN_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Min : protected Pointers {
 public:
  double einitial, efinal, eprevious;
  double fnorm2_init, fnorminf_init, fnorm2_final, fnorminf_final;
  double alpha_final;
  int niter, neval;
  int stop_condition;
  char *stopstr;
  int searchflag;    // 0 if damped dynamics, 1 if sub-cycles on local search

  Min(class LAMMPS *);
  ~Min() override;
  virtual void init();
  virtual void setup(int flag = 1);
  virtual void setup_minimal(int);
  virtual void run(int);
  virtual void force_clear();
  void cleanup();
  int request(class Pair *, int, double);
  virtual double memory_usage() { return 0; }
  void modify_params(int, char **);
  virtual int modify_param(int, char **) { return 0; }
  virtual double fnorm_sqr();
  virtual double fnorm_inf();
  virtual double fnorm_max();

  enum { TWO, MAX, INF };

  // methods for spin minimizers
  double total_torque();
  double inf_torque();
  double max_torque();

  virtual void init_style() {}
  virtual void setup_style() = 0;
  virtual void reset_vectors() = 0;
  virtual int iterate(int) = 0;

  // possible return values of iterate() method
  enum {
    MAXITER,
    MAXEVAL,
    ETOL,
    FTOL,
    DOWNHILL,
    ZEROALPHA,
    ZEROFORCE,
    ZEROQUAD,
    TRSMALL,
    INTERROR,
    TIMEOUT,
    MAXVDOTF
  };

  // integrator styles
  enum { EULERIMPLICIT, VERLET, LEAPFROG, EULEREXPLICIT };

  // line search styles
  enum { BACKTRACK, QUADRATIC, FORCEZERO, SPIN_CUBIC, SPIN_NONE };

 protected:
  int eflag, vflag;            // flags for energy/virial computation
  int virial_style;            // compute virial explicitly or implicitly
  int external_force_clear;    // clear forces locally or externally

  double dmax;      // max dist to move any atom in one step
  int linestyle;    // 0 = backtrack, 1 = quadratic, 2 = forcezero
                    // 3 = spin_cubic, 4 = spin_none

  int normstyle;    // TWO, MAX or INF flag for force norm evaluation

  double dtinit;    // store the default timestep

  // only for minimize style fire2
  int delaystep;                 // minium steps of dynamics
  double dtgrow, dtshrink;       // timestep increase, decrease
  double alpha0, alphashrink;    // mixing velocities+forces coefficient
  double tmax, tmin;             // timestep multiplicators max, min
  int integrator;                // Newton integration: euler, leapfrog, verlet...
  int halfstepback_flag;         // half step backward when v.f <= 0.0
  int delaystep_start_flag;      // delay the initial dt_shrink
  int max_vdotf_negatif;         // maximum iteration with v.f > 0.0
  int abcflag;                   // when 1 use ABC-FIRE variant instead of FIRE, default 0

  int nelist_global, nelist_atom;    // # of PE,virial computes to check
  int nvlist_global, nvlist_atom, ncvlist_atom;
  class Compute **elist_global;    // lists of PE,virial Computes
  class Compute **elist_atom;
  class Compute **vlist_global;
  class Compute **vlist_atom;
  class Compute **cvlist_atom;

  int triclinic;    // 0 if domain is orthog, 1 if triclinic
  int pairflag;
  int torqueflag, extraflag;

  int pair_compute_flag;      // 0 if pair->compute is skipped
  int kspace_compute_flag;    // 0 if kspace->compute is skipped

  int narray;                         // # of arrays stored by fix_minimize
  class FixMinimize *fix_minimize;    // fix that stores auxiliary data

  class Compute *pe_compute;    // compute for potential energy
  double ecurrent;              // current potential energy

  bigint ndoftotal;    // total dof for entire problem

  int nvec;        // local atomic dof = length of xvec
  double *xvec;    // variables for atomic dof, as 1d vector
  double *fvec;    // force vector for atomic dof, as 1d vector

  int nextra_global;    // # of extra global dof due to fixes
  double *fextra;       // force vector for extra global dof
                        // xextra is stored by fix

  int nextra_atom;           // # of extra per-atom variables
  double **xextra_atom;      // ptr to the variable
  double **fextra_atom;      // ptr to the force on the variable
  int *extra_peratom;        // # of values in variable, e.g. 3 in x
  int *extra_nlen;           // total local length of variable, e.g 3*nlocal
  double *extra_max;         // max allowed change per iter for atom's var
  class Pair **requestor;    // Pair that stores/manipulates the variable

  int kokkosable;    // 1 if this min style supports Kokkos

  int neigh_every, neigh_delay, neigh_dist_check;    // neighboring params

  virtual double energy_force(int);

  void ev_setup();
  void ev_set(bigint);

  char *stopstrings(int);
};

}    // namespace LAMMPS_NS

#endif
