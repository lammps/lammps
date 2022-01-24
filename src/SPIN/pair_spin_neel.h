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
PairStyle(spin/neel,PairSpinNeel);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_NEEL_H
#define LMP_PAIR_SPIN_NEEL_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinNeel : public PairSpin {
 public:
  PairSpinNeel(LAMMPS *lmp) : PairSpin(lmp) {}
  virtual ~PairSpinNeel();
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_neel(int, int, double, double *, double *, double *, double *);
  void compute_neel_mech(int, int, double, double *, double *, double *, double *);
  double compute_neel_energy(int, int, double, double *, double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  double cut_spin_neel_global;    // global neel cutoff distance

 protected:
  // pseudo-dipolar and pseudo-quadrupolar coeff.

  double **g1, **g1_mech;    // neel coeffs gij
  double **g2, **g3;         // g1 in eV, g2 adim, g3 in Ang
  double **q1, **q1_mech;    // neel coeffs qij
  double **q2, **q3;         // q1 in eV, q2 adim, q3 in Ang
  double **cut_spin_neel;    // cutoff distance exchange

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_spin command

Self-explanatory.

E: Spin simulations require metal unit style

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair spin requires atom attribute spin

The atom style defined does not have these attributes.

*/
