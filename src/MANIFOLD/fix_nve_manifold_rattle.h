/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterized by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/manifold/rattle,FixNVEManifoldRattle);
// clang-format on
#else

#ifndef LMP_FIX_NVE_MANIFOLD_RATTLE_H
#define LMP_FIX_NVE_MANIFOLD_RATTLE_H

#include "fix.h"

namespace LAMMPS_NS {

namespace user_manifold {
  class manifold;
}

class FixNVEManifoldRattle : public Fix {
 public:
  struct statistics {

    statistics() :
        x_iters(0), v_iters(0), x_iters_per_atom(0), v_iters_per_atom(0), natoms(0),
        dofs_removed(0), last_out(0)
    {
    }
    double x_iters, v_iters;
    double x_iters_per_atom;
    double v_iters_per_atom;
    int natoms;
    int dofs_removed;
    bigint last_out;
  };

  FixNVEManifoldRattle(LAMMPS *, int &, char **, int = 1);
  ~FixNVEManifoldRattle() override;
  // All this stuff is interface, so you DO need to implement them.
  // Just delegate them to the workhorse classes.
  int setmask() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void init() override;
  void reset_dt() override;
  void end_of_step() override;
  bigint dof(int) override;
  void setup(int) override {}    // Not needed for fixNVE but is for fixNVT
  double memory_usage() override;

 protected:
  int nevery, next_output;

  double dtv, dtf;
  double tolerance;
  int max_iter;

  char **tstrs;
  int *tvars;
  int *tstyle;
  int *is_var;

  statistics stats;
  int update_style;
  int nvars;

  user_manifold::manifold *ptr_m;

  void print_stats(const char *);
  int was_var(const char *);

  virtual void update_var_params();
  virtual void rattle_manifold_x(double *, double *, double *, double, double, tagint);
  virtual void rattle_manifold_v(double *, double *, double *, double, tagint);

  virtual void nve_x_rattle(int, int);
  virtual void nve_v_rattle(int, int);
};
}    // namespace LAMMPS_NS

#endif    // LMP_FIX_NVE_MANIFOLD_RATTLE_H
#endif
