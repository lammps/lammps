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
PairStyle(lj/long/dipole/long,PairLJLongDipoleLong);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_LONG_DIPOLE_LONG_H
#define LMP_PAIR_LJ_LONG_DIPOLE_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJLongDipoleLong : public Pair {
 public:
  double cut_coul;

  PairLJLongDipoleLong(class LAMMPS *);
  virtual ~PairLJLongDipoleLong();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);

  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void *extract(const char *, int &);

 protected:
  double cut_lj_global;
  double **cut_lj, **cut_lj_read, **cut_ljsq;
  double cut_coulsq;
  double **epsilon_read, **epsilon, **sigma_read, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;
  double *cut_respa;
  double g_ewald;
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

W: Geometric mixing assumed for 1/r^6 coefficients

Self-explanatory.

W: Using largest cut-off for lj/long/dipole/long long long

Self-explanatory.

E: Cut-offs missing in pair_style lj/long/dipole/long

Self-explanatory.

E: Coulombic cut not supported in pair_style lj/long/dipole/long

Must use long-range Coulombic interactions.

E: Only one cut-off allowed when requesting all long

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Invoking coulombic in pair style lj/long/dipole/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair lj/long/dipole/long requires atom attributes mu, torque

The atom style defined does not have these attributes.

E: Pair style requires a KSpace style

No kspace style is defined.

E: Pair style requires use of kspace_style ewald/disp

Self-explanatory.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

U: Pair style lj/long/dipole/long does not currently support respa

This feature is not yet supported.

*/
