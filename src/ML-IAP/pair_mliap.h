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
PairStyle(mliap,PairMLIAP);
// clang-format on
#else

#ifndef LMP_PAIR_MLIAP_H
#define LMP_PAIR_MLIAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMLIAP : public Pair {
 public:
  PairMLIAP(class LAMMPS *);
  ~PairMLIAP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void e_tally(class MLIAPData *);
  void v_tally(int, int, double *, double *);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual double memory_usage();
  int *map;    // mapping from atom types to elements

 protected:
  virtual void allocate();

  class MLIAPModel *model;
  class MLIAPDescriptor *descriptor;
  class MLIAPData *data;
};

}    // namespace LAMMPS_NS

#endif
#endif
