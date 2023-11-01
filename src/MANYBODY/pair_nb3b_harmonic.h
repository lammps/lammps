/* -*- c++ -*- ---------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(nb3b/harmonic,PairNb3bHarmonic);
// clang-format on
#else

#ifndef LMP_PAIR_NB3B_HARMONIC_H
#define LMP_PAIR_NB3B_HARMONIC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNb3bHarmonic : public Pair {
 public:
  PairNb3bHarmonic(class LAMMPS *);
  ~PairNb3bHarmonic() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  static constexpr int NPARAMS_PER_LINE = 6;
  enum { HARMONIC = 0, SCREENED };

 protected:
  struct Param {
    double k_theta, theta0, cutoff;
    double invrho;    // for screened harmonic style
    double cut, cutsq;
    int ielement, jelement, kelement;
  };

  double cutmax;    // max cutoff for all elements
  Param *params;    // parameter set for an I-J-K interaction
  int variant;

  void allocate();
  void read_file(char *);
  void setup_params();
  void twobody(Param *, double, double &, int, double &);
  virtual void threebody(Param *, Param *, Param *, double, double, double *, double *, double *,
                         double *, int, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif
