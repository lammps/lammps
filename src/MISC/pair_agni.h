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
PairStyle(agni,PairAGNI);
// clang-format on
#else

#ifndef LMP_PAIR_AGNI_H
#define LMP_PAIR_AGNI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairAGNI : public Pair {
 public:
  enum { AGNI_VERSION_UNKNOWN, AGNI_VERSION_1, AGNI_VERSION_2 };
  PairAGNI(class LAMMPS *);
  ~PairAGNI() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  struct Param {
    double cut, cutsq;
    double *eta, **xU, *alpha;
    double sigma, lambda, b, gwidth;
    int numeta, numtrain, ielement;
  };

 protected:
  double cutmax;                 // max cutoff for all elements
  int atomic_feature_version;    // version of fingerprint
  Param *params;                 // parameter set for an I-J interaction
  virtual void allocate();
  void read_file(char *);
  virtual void setup_params();
};

}    // namespace LAMMPS_NS

#endif
#endif
