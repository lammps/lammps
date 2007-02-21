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
  void init();
  void run();

  virtual void iterate(int);

 protected:
  int virial_thermo;         // what vflag should be on thermo steps (1,2)
  int pairflag,torqueflag,granflag;               // force clear flags
  int neigh_every,neigh_delay,neigh_dist_check;   // copies of reneigh criteria

  int maxpair;               // copies of Update quantities
  double **f_pair;

  class FixMinimize *fix_minimize;  // fix that stores gradient vecs
  double mindist,maxdist;    // min/max dist for coord delta in line search

  int ndof;                  // 3N degrees-of-freedom on this proc
  double *x;                 // vec of 3N coords, ptr to atom->x[0]
  double *f;                 // vec of 3N forces, ptr to atom->f[0]
  double *g;                 // vec of 3N old forces, ptr to fix_minimize::g
  double *h;                 // vec of 3N search dir, ptr to fix_minimize::h

  int ndof_extra;            // extra degrees of freedom to include in min
  double energy_extra;       // extra potential energy
  double *xextra;            // extra vecs associated with ndof_extra
  double *fextra;
  double *gextra;
  double *hextra;

  double energy;             // potential energy of atoms and extra dof

  typedef int (MinCG::*FnPtr)(int &);
  FnPtr linemin;                        // ptr to linemin functions
  int linemin_scan(int &);
  int linemin_secant(int &);

  void setup();
  void set_local_vectors();
  void eng_force();
  void force_clear(int);
};

}

#endif
