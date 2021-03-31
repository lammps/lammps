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

#ifdef MINIMIZE_CLASS

MinimizeStyle(spin/cg, MinSpinCG)

#else

#ifndef LMP_MIN_SPIN_CG_H
#define LMP_MIN_SPIN_CG_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpinCG: public Min {
 public:
  MinSpinCG(class LAMMPS *);
  virtual ~MinSpinCG();
  void init();
  void setup_style();
  void reset_vectors();
  int modify_param(int, char **);
  int iterate(int);

 private:
  int local_iter;               // for neb
  int nlocal_max;               // max value of nlocal (for size of lists)
  int use_line_search;          // use line search or not.
  int ireplica,nreplica;        // for neb
  double dt;                    // global timestep
  double dts;                   // spin timestep
  double discrete_factor;       // factor for spin timestep evaluation
  double der_e_cur;             // current derivative along search dir.
  double der_e_pr;              // previous derivative along search dir.
  double *spvec;                // variables for atomic dof, as 1d vector
  double *fmvec;                // variables for atomic dof, as 1d vector
  double *g_old;                // gradient vector at previous step
  double *g_cur;                // current gradient vector
  double *p_s;                  // search direction vector
  double **sp_copy;             // copy of the spins

  void advance_spins();
  void calc_gradient();
  void calc_search_direction();
  void vm3(const double *, const double *, double *);
  void rodrigues_rotation(const double *, double *);
  void make_step(double, double *);
  int calc_and_make_step(double, double, int);
  int adescent(double, double);
  double evaluate_dt();
  double maximum_rotation(double *);

  bigint last_negative;
};

}

#endif
#endif
