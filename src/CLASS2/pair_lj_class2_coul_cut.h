/* ----------------------------------------------------------------------
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

PairStyle(lj/class2/coul/cut,PairLJClass2CoulCut)

#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_CUT_H
#define LMP_PAIR_LJ_CLASS2_COUL_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJClass2CoulCut : public Pair {
 public:
  PairLJClass2CoulCut(class LAMMPS *);
  virtual ~PairLJClass2CoulCut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_lj_global,cut_coul_global;
  double **cut_lj,**cut_ljsq;
  double **cut_coul,**cut_coulsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

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

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/class2/coul/cut requires atom attribute q

The atom style defined does not have this attribute.

*/
