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

MinimizeStyle(spin_oso_lbfgs, MinSpinOSO_LBFGS)

#else

#ifndef LMP_MIN_SPIN_OSO_LBFGS_H
#define LMP_MIN_SPIN_OSO_LBFGS_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpinOSO_LBFGS: public Min {
 public:
  MinSpinOSO_LBFGS(class LAMMPS *);
  virtual ~MinSpinOSO_LBFGS();
  void init();
  void setup_style();
  int modify_param(int, char **);
  void reset_vectors();
  int iterate(int);

 private:
  int local_iter;     		// for neb
  int use_line_search; 		// use line search or not.
  int nlocal_max;		// max value of nlocal (for size of lists)
  int ireplica,nreplica; 	// for neb
  double der_e_cur;		// current derivative along search dir.
  double der_e_pr;    		// previous derivative along search dir.
  double maxepsrot;
  double *spvec;		// variables for atomic dof, as 1d vector
  double *fmvec;		// variables for atomic dof, as 1d vector
  double *g_old;  		// gradient vector at previous step
  double *g_cur;  		// current gradient vector
  double *p_s;  		// search direction vector

  void advance_spins();
  void calc_gradient();
  void calc_search_direction();
  void vm3(const double *, const double *, double *);
  void rodrigues_rotation(const double *, double *);
  void make_step(double, double *);
  int calc_and_make_step(double, double, int);
  int adescent(double, double);
  double maximum_rotation(double *);

  double *rho;		// estimation of curvature
  double **ds;		// change in rotation matrix between two iterations, da
  double **dy;		// change in gradients between two iterations, dg
  double **sp_copy;	// copy of the spins
  int num_mem;		// number of stored steps
  bigint last_negative;
};

}

#endif
#endif
