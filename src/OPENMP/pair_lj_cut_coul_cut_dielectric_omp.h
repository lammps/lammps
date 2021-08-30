/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/cut/dielectric/omp,PairLJCutCoulCutDielectricOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_DIELECTRIC_OMP_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_DIELECTRIC_OMP_H

#include "pair_lj_cut_coul_cut_dielectric.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCutCoulCutDielectricOMP : public PairLJCutCoulCutDielectric, public ThrOMP {
 public:
  PairLJCutCoulCutDielectricOMP(class LAMMPS *);
  virtual ~PairLJCutCoulCutDielectricOMP();
  virtual void compute(int, int);

 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/cut/coul/cut/dielectric requires atom attribute q

The atom style defined does not have this attribute.

*/
