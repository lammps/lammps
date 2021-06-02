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

#ifdef FIX_CLASS
// clang-format off
FixStyle(latte,FixLatte);
// clang-format on
#else

#ifndef LMP_FIX_LATTE_H
#define LMP_FIX_LATTE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLatte : public Fix {
 public:
  FixLatte(class LAMMPS *, int, char **);
  virtual ~FixLatte();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void min_setup(int);
  void setup_pre_reverse(int, int);
  void initial_integrate(int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  void final_integrate();
  void reset_dt();
  double compute_scalar();
  double memory_usage();

 protected:
  char *id_pe;
  int coulomb, pbcflag, pe_peratom, virial_global, virial_peratom, neighflag;
  int eflag_caller;

  int nmax, newsystem;
  double *qpotential;
  double **flatte;
  double latte_energy;

  class NeighList *list;
  class Compute *c_pe;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must use units metal with fix latte command

UNDOCUMENTED

E: Fix latte currently runs only in serial

UNDOCUMENTED

E: LAMMPS is linked against incompatible LATTE library

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix latte does not yet support a LAMMPS calculation of a Coulomb potential

UNDOCUMENTED

E: Could not find fix latte compute ID

UNDOCUMENTED

E: Fix latte compute ID does not compute pe/atom

UNDOCUMENTED

E: Fix latte requires 3d problem

UNDOCUMENTED

E: Fix latte cannot compute Coulomb potential

UNDOCUMENTED

E: Fix latte requires 3d simulation

UNDOCUMENTED

W: Fix latte should come after all other integration fixes

UNDOCUMENTED

E: Internal LATTE problem

UNDOCUMENTED

*/
