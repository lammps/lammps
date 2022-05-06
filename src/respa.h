/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS
// clang-format off
IntegrateStyle(respa,Respa);
// clang-format on
#else

#ifndef LMP_RESPA_H
#define LMP_RESPA_H

#include "integrate.h"

namespace LAMMPS_NS {

class Respa : public Integrate {
 public:
  // public so Fixes, Pairs, Neighbor can see them
  int nlevels;         // number of rRESPA levels
                       // 0 = innermost level, nlevels-1 = outermost level
  double *step;        // timestep at each level
  int *loop;           // sub-cycling factor at each level
  double cutoff[4];    // cutoff[0] and cutoff[1] = between inner and middle
                       // cutoff[2] and cutoff[3] = between middle and outer
                       // if no middle then 0,1 = 2,3

  int level_bond, level_angle, level_dihedral;    // level to compute forces at
  int level_improper, level_pair, level_kspace;
  int level_inner, level_middle, level_outer;

  int nhybrid_styles;     // number of hybrid pair styles
  int *hybrid_level;      // level to compute pair hybrid sub-style at
  int *hybrid_compute;    // selects whether to compute sub-style forces
  int tally_global;       // 1 if pair style should tally global accumulators
  int pair_compute;       // 1 if pair force need to be computed

  Respa(class LAMMPS *, int, char **);
  ~Respa() override;
  void init() override;
  void setup(int) override;
  void setup_minimal(int) override;
  void run(int) override;
  void force_clear() override;
  void cleanup() override;
  void reset_dt() override;

  void copy_f_flevel(int);
  void copy_flevel_f(int);

 protected:
  int triclinic;    // 0 if domain is orthog, 1 if triclinic
  int torqueflag, extraflag;

  int *newton;                  // newton flag at each level
  class FixRespa *fix_respa;    // Fix to store the force level array

  virtual void recurse(int);
  void sum_flevel_f();
  void set_compute_flags(int ilevel);
};

}    // namespace LAMMPS_NS

#endif
#endif
