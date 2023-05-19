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
// clang-format off
FixStyle(adapt/dcci,FixAdaptDCCI);
// clang-format on
#else

#ifndef LMP_FIX_ADAPT_DCCI_H
#define LMP_FIX_ADAPT_DCCI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdaptDCCI : public Fix {
 public:
  FixAdaptDCCI(class LAMMPS *, int, char **);
  ~FixAdaptDCCI() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void post_run() override;
  void setup_pre_force_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  double lambda;
  double compute_scalar() override;

 private:
  int nadapt, fscaleflag;
  int anypair;
  int nlevels_respa;

  struct Adapt {
    int which;
    char *pstyle, *pparam;
    int ilo, ihi, jlo, jhi;
    int pdim;
    double *scalar, scalar_orig;
    double **array, **array_orig;
    class Pair *pair;
  };

  Adapt *adapt;

  void change_settings();
  void restore_settings();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix adapt/dcci does not exist

Self-explanatory.

E: Variable for fix adapt/dcci is invalid style

Only equal-style variables can be used.

E: Fix adapt/dcci pair style does not exist

Self-explanatory

E: Fix adapt/dcci pair style param not supported

The pair style does not know about the parameter you specified.

E: Fix adapt/dcci pair style param is not compatible

Self-explanatory

E: Fix adapt/dcci type pair range is not valid for pair hybrid sub-style

Self-explanatory.

E: Fix adapt/dcci kspace style does not exist

Self-explanatory.

E: Could not find fix adapt/dcci storage fix ID

This should not happen unless you explicitly deleted
a secondary fix that fix adapt created internally.

*/
