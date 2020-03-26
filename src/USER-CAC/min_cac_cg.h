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

MinimizeStyle(cac/cg,CACMinCG)

#else

#ifndef LMP_CAC_MIN_CG_H
#define LMP_CAC_MIN_CG_H

#include "min.h"

namespace LAMMPS_NS {

class CACMinCG : public Min {
 public:
	 CACMinCG(class LAMMPS *);
  ~CACMinCG();
  void init();
  void setup_style();
  void reset_vectors();
  virtual int iterate(int);
 protected:
  // vectors needed by linesearch minimizers
  // allocated and stored by fix_minimize
  // x,f are stored by parent or Atom class or Pair class

  double *x0;                 // coords at start of linesearch
  double *g;                  // old gradient vector
  double *h;                  // search direction vector

  double *gextra;             // g,h for extra global dof, x0 is stored by fix
  double *hextra;

  double **x0extra_atom;      // x0,g,h for extra per-atom dof
  double **gextra_atom;
  double **hextra_atom;

  int densemax;               // bounds arrays size for continuous x,v,f nodal arrays


  typedef int (CACMinCG::*FnPtr)(double, double &);
  FnPtr linemin;
  int linemin_backtrack(double, double &);
  int linemin_quadratic(double, double &);
  int linemin_forcezero(double, double &);

  double alpha_step(double, int);
  double compute_dir_deriv(double &);
  virtual void copy_vectors();
  virtual void copy_force();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: CAC min styles require a CAC atom style

Self Explanatory

E: CAC min styles require a CAC pair style

Self Explanatory

*/