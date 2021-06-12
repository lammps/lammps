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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(buck/long/coul/long/omp,PairBuckLongCoulLongOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BUCK_LONG_COUL_LONG_OMP_H
#define LMP_PAIR_BUCK_LONG_COUL_LONG_OMP_H

#include "pair_buck_long_coul_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairBuckLongCoulLongOMP : public PairBuckLongCoulLong, public ThrOMP {

 public:
  PairBuckLongCoulLongOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual void compute_inner();
  virtual void compute_middle();
  virtual void compute_outer(int, int);

 private:
  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval(int, int, ThrData *const);

  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval_outer(int, int, ThrData *const);

  void eval_inner(int, int, ThrData *const);
  void eval_middle(int, int, ThrData *const);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Geometric mixing assumed for 1/r^6 coefficients

Self-explanatory.

W: Using largest cutoff for buck/long/coul/long

Self-explanatory.

E: Cutoffs missing in pair_style buck/long/coul/long

Self-explanatory.

E: LJ6 off not supported in pair_style buck/long/coul/long

Self-explanatory.

E: Coulomb cut not supported in pair_style buck/long/coul/coul

Must use long-range Coulombic interactions.

E: Only one cutoff allowed when requesting all long

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style buck/long/coul/long requires atom attribute q

The atom style defined does not have this attribute.

E: Pair style requires a KSpace style

No kspace style is defined.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
