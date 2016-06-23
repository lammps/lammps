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
/* ----------------------------------------------------------------------
   Contributing author: Stefan Paquay (Eindhoven University of Technology)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(morse/soft,PairMorseSoft)

#else

#ifndef LMP_PAIR_MORSE_SOFT_H
#define LMP_PAIR_MORSE_SOFT_H

#include "pair_morse.h"

namespace LAMMPS_NS {

class PairMorseSoft : public PairMorse {
 public:
  PairMorseSoft(class LAMMPS *lmp) : PairMorse(lmp) {};
  virtual ~PairMorseSoft();
  virtual void compute(int, int);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);

  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double **lambda;

  int nlambda;
  double shift_range;

  void allocate();
};

}

#endif
#endif
