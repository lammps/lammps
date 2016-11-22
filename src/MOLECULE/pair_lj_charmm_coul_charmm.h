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

PairStyle(lj/charmm/coul/charmm,PairLJCharmmCoulCharmm)

#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_CHARMM_H
#define LMP_PAIR_LJ_CHARMM_COUL_CHARMM_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulCharmm : public Pair {
 public:
  PairLJCharmmCoulCharmm(class LAMMPS *);
  virtual ~PairLJCharmmCoulCharmm();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

 protected:
  int implicit;
  double cut_lj_inner,cut_lj,cut_coul_inner,cut_coul;
  double cut_lj_innersq,cut_ljsq,cut_coul_innersq,cut_coulsq,cut_bothsq;
  double denom_lj,denom_coul;
  double **epsilon,**sigma,**eps14,**sigma14;
  double **lj1,**lj2,**lj3,**lj4;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/charmm/coul/charmm requires atom attribute q

The atom style defined does not have these attributes.

E: Pair inner cutoff >= Pair outer cutoff

The specified cutoffs for the pair style are inconsistent.

*/
