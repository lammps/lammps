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

#ifdef FIX_CLASS

FixStyle(ti/rs,FixTIRS)

#else

#ifndef LMP_FIX_TI_RS_H
#define LMP_FIX_TI_RS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTIRS : public Fix {
 public:
  FixTIRS(class LAMMPS *, int, char **);
  ~FixTIRS();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  virtual void initial_integrate(int);
  double compute_vector(int);

 private:
  double switch_func(double);
  double dswitch_func(double);

  double lambda;       // Coupling parameter.
  double dlambda;      // Lambda variation with t.
  double l_initial;    // Lambda initial value.
  double l_final;      // Lambda final value.
  double linfo[2];     // Current lambda status.
  int    t_switch;     // Total switching steps.
  int    t_equil;      // Equilibration time.
  int    t0;           // Initial time.
  int    sf;           // Switching function option.
  int    nlevels_respa;
};

}

#endif
#endif
