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
PairStyle(lj/long/coul/long,PairLJLongCoulLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJLongCoulLong : public Pair {
 public:
  double cut_coul;

  PairLJLongCoulLong(class LAMMPS *);
  ~PairLJLongCoulLong() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;

 protected:
  double cut_lj_global;
  double **cut_lj, **cut_lj_read, **cut_ljsq;
  double cut_coulsq;
  double **epsilon_read, **epsilon, **sigma_read, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;
  double qdist;
  double g_ewald;
  double g_ewald_6;
  int ewald_order, ewald_off;

  void options(char **arg, int order);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Using largest cutoff for lj/long/coul/long

Self-explanatory.

E: Cutoffs missing in pair_style lj/long/coul/long

Self-explanatory.

E: Coulomb cut not supported in pair_style lj/long/coul/long

Must use long-range Coulombic interactions.

E: Only one cutoff allowed when requesting all long

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Invoking coulombic in pair style lj/long/coul/long requires atom attribute q

UNDOCUMENTED

E: Pair style requires a KSpace style

No kspace style is defined.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

U: Invoking coulombic in pair style lj/coul requires atom attribute q

The atom style defined does not have this attribute.

*/
