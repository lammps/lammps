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

#ifdef FIX_CLASS

FixStyle(grem,FixGrem)

#else

#ifndef LMP_FIX_GREM_H
#define LMP_FIX_GREM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGrem : public Fix {
 public:
  FixGrem(class LAMMPS *, int, char **);
  ~FixGrem();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void *extract(const char *, int &);
  double compute_scalar();
  double scale_grem,lambda,eta,h0;
  int pressflag;

 private:
  double tbath,pressref;

 protected:
  char *id_temp,*id_press,*id_ke,*id_pe,*id_nh;
  class Compute *temperature,*pressure,*ke,*pe;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix grem does not exist

Self-explanatory.

E: Variable name for fix grem does not exist

Self-explanatory.

E: Variable for fix grem is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix grem

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix grem

Must define an energy variable when applying a dynamic
force during minimization.

*/
