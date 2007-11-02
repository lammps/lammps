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

#ifndef MIN_CG_H
#define MIN_CG_H

#include "min.h"

namespace LAMMPS_NS {

class MinCG : public Min {
 public:
  MinCG(class LAMMPS *);
  virtual ~MinCG();
  void init();
  void run();
  virtual void iterate(int);

 protected:
  int eflag,vflag;            // flags for energy/virial computation
  int virial_style;           // compute virial explicitly or implicitly
  int pairflag,torqueflag;
  int neigh_every,neigh_delay,neigh_dist_check;   // copies of reneigh criteria
  int triclinic;              // 0 if domain is orthog, 1 if triclinic

  int nvlist;                 // # of PE,virial coputes for eflag,vflag
  class Compute **vlist;      // list of Computes to check

  int maxpair;                // copies of Update quantities
  double **f_pair;

  class FixMinimize *fix_minimize;  // fix that stores gradient vecs
  class Compute *pe_compute;  // compute for potential energy
  double ecurrent;            // current potential energy
  double mindist,maxdist;     // min/max dist for coord delta in line search

  int ndof;                   // # of degrees-of-freedom on this proc
  double *g,*h;               // local portion of gradient, searchdir vectors

  typedef int (MinCG::*FnPtr)(int, double *, double *, double,
			      double, double, double &, int &);
  FnPtr linemin;             // ptr to linemin functions

  int linemin_scan(int, double *, double *, double,
		   double, double, double &, int &);
  int linemin_secant(int, double *, double *, double,
		     double, double, double &, int &);

  void setup();
  void setup_vectors();
  void eng_force(int *, double **, double **, double *);
  void ev_set(int);
  void force_clear(int);
};

}

#endif
