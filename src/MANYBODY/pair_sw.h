/* -*- c++ -*- ----------------------------------------------------------
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
PairStyle(sw,PairSW);
// clang-format on
#else

#ifndef LMP_PAIR_SW_H
#define LMP_PAIR_SW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSW : public Pair {
 public:
  PairSW(class LAMMPS *);
  ~PairSW() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  double single(int, int, int, int, double, double, double, double &) override;

  static constexpr int NPARAMS_PER_LINE = 14;

  struct Param {
    double epsilon, sigma;
    double littlea, lambda, gamma, costheta;
    double biga, bigb;
    double powerp, powerq;
    double tol;
    double cut, cutsq;
    double sigma_gamma, lambda_epsilon, lambda_epsilon2;
    double c1, c2, c3, c4, c5, c6;
    int ielement, jelement, kelement;
  };

 protected:
  double cutmax;              // max cutoff for all elements
  Param *params;              // parameter set for an I-J-K interaction
  int maxshort;               // size of short neighbor list array
  int *neighshort;            // short neighbor list array
  int skip_threebody_flag;    // whether to run threebody loop
  int params_mapped;          // whether parameters have been read and mapped to elements

  void settings(int, char **) override;
  virtual void allocate();
  virtual void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);
  virtual void threebody(Param *, Param *, Param *, double, double, double *, double *, double *,
                         double *, int, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif
