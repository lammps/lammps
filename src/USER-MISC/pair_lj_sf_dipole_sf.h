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
PairStyle(lj/sf/dipole/sf,PairLJSFDipoleSF);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_SF_DIPOLE_SF_H
#define LMP_PAIR_LJ_SF_DIPOLE_SF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSFDipoleSF : public Pair {
 public:
  PairLJSFDipoleSF(class LAMMPS *);
  virtual ~PairLJSFDipoleSF();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_lj_global, cut_coul_global;
  double **cut_lj, **cut_ljsq;
  double **cut_coul, **cut_coulsq;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4;
  double **scale;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
