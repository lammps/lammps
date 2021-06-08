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
  virtual ~PairAGNI();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open AGNI potential file %s

The specified AGNI potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in AGNI potential file

The potential file is not compatible with the AGNI pair style
implementation in this LAMMPS version.

*/
