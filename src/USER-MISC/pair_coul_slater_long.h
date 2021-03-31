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

PairStyle(coul/slater/long,PairCoulSlaterLong)

#else

#ifndef LMP_PAIR_COUL_SLATER_LONG_H
#define LMP_PAIR_COUL_SLATER_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulSlaterLong : public Pair {
 public:
  PairCoulSlaterLong(class LAMMPS *);
  ~PairCoulSlaterLong();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

 protected:
  double cut_coul,cut_coulsq,qdist;
  double lamda;
  //double *cut_respa;
  double g_ewald;
  double **scale;

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

E: Pair style coul/slater/long requires atom attribute q

The atom style defined does not have this attribute.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
