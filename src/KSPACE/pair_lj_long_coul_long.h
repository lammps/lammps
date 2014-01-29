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

#ifdef PAIR_CLASS

PairStyle(lj/long/coul/long,PairLJLongCoulLong)

#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_H
 
#include "pair.h"

namespace LAMMPS_NS {

class PairLJLongCoulLong : public Pair {
 public:
  double cut_coul;

  PairLJLongCoulLong(class LAMMPS *);
  virtual ~PairLJLongCoulLong();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  virtual void compute_inner();
  virtual void compute_middle();
  virtual void compute_outer(int, int);

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

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Using largest cutoff for lj/long/coul/long

UNDOCUMENTED

E: Cutoffs missing in pair_style lj/long/coul/long

Self-explanatory.

E: Coulomb cut not supported in pair_style lj/long/coul/long

Must use long-range Coulombic interactions.

E: Only one cutoff allowed when requesting all long

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Invoking coulombic in pair style lj/coul requires atom attribute q

UNDOCUMENTED

E: Pair style requires a KSpace style

No kspace style is defined.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

U: Mixing forced for lj coefficients

Self-explanatory.

U: Mixing forced for LJ coefficients

Self-explanatory.

U: Using largest cutoff for pair_style lj/long/coul/long

Self-explanatory.

U: Pair style lj/long/coul/long requires atom attribute q

The atom style defined does not have this attribute.

*/
