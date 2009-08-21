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

#ifndef RESPA_H
#define RESPA_H

#include "integrate.h"

namespace LAMMPS_NS {

class Respa : public Integrate {
 public:
                          // public so Fixes, Pairs, Neighbor can see them
  int nlevels;            // number of rRESPA levels
                          // 0 = innermost level, nlevels-1 = outermost level
  double *step;           // timestep at each level
  int *loop;              // sub-cycling factor at each level
  double cutoff[4];       // cutoff[0] and cutoff[1] = between inner and middle
                          // cutoff[2] and cutoff[3] = between middle and outer
                          // if no middle then 0,1 = 2,3

  int level_bond,level_angle,level_dihedral;   // level to compute forces at
  int level_improper,level_pair,level_kspace;
  int level_inner,level_middle,level_outer;

  Respa(class LAMMPS *, int, char **);
  ~Respa();
  void init();
  void setup();
  void setup_minimal(int);
  void run(int);
  void cleanup();
  void reset_dt();

  void copy_f_flevel(int);
  void copy_flevel_f(int);

 private:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic

  int *newton;                      // newton flag at each level
  class FixRespa *fix_respa;        // Fix to store the force level array

  void recurse(int);
  void force_clear(int);
  void sum_flevel_f();
};

}

#endif
