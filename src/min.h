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
  void run();
  double memory_usage() {return 0.0;}
  void modify_params(int, char **);

  virtual int iterate(int) = 0;

 protected:
  int eflag,vflag;            // flags for energy/virial computation
  int virial_style;           // compute virial explicitly or implicitly

  double dmax;                // max dist to move any atom in one linesearch
  int linestyle;              // 0 = backtrack, 1 = quadratic

  int nelist_atom;                  // # of PE,virial computes to check
  int nvlist_global,nvlist_atom;
  class Compute **elist_atom;       // list of PE,virial Computes
  class Compute **vlist_global;
  class Compute **vlist_atom;

  int pairflag,torqueflag;
  int neigh_every,neigh_delay,neigh_dist_check;   // copies of reneigh criteria
  int triclinic;              // 0 if domain is orthog, 1 if triclinic

  class FixMinimize *fix_minimize;  // fix that stores gradient vecs
  class Compute *pe_compute;        // compute for potential energy
  double ecurrent;                  // current potential energy
  double mindist,maxdist;     // min/max dist for coord delta in line search

  int ndof;                   // # of degrees-of-freedom on this proc
  double *g,*h;               // local portion of gradient, searchdir vectors
  double *x0;                 // coords at start of linesearch

  int nextra;                 // extra dof due to fixes
  double *fextra;             // vectors for extra dof
  double *gextra;
  double *hextra;

  double boxlo0[3];           // box size at start of linesearch
  double boxhi0[3];

  // ptr to linemin functions

  void setup();
  void eng_force(int *, double **, double **, double **, double *, int);
  void setup_vectors();
  void force_clear();

  typedef int (Min::*FnPtr)(int, double *, double *, double *, double,
			    double, double &, int &);
  FnPtr linemin;
  int linemin_backtrack(int, double *, double *, double *, double,
			double, double &, int &);
  int linemin_quadratic(int, double *, double *, double *, double,
			double, double &, int &);

  void ev_setup();
  void ev_set(int);
  void box_store();
  void box_swap();
};

}

#endif
