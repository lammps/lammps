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
PairStyle(edip/multi,PairEDIPMulti);
// clang-format on
#else

#ifndef LMP_PAIR_EDIP_MULTI_H
#define LMP_PAIR_EDIP_MULTI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEDIPMulti : public Pair {
 public:
  PairEDIPMulti(class LAMMPS *);
  ~PairEDIPMulti() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  static constexpr int NPARAMS_PER_LINE = 20;

 protected:
  struct Param {
    double A, B;              // coefficients for pair interaction I-J
    double cutoffA;           // cut-off distance for pair interaction I-J
    double cutoffC;           // lower cut-off distance for calculating Z_I
    double alpha;             // coefficient for calculating Z_I
    double beta;              // attractive term for pair I-J
    double sigma;             // cut-off coefficient for pair I-J
    double rho;               // pair I-J
    double gamma;             // coefficient for three-body interaction I-J-K
    double eta, lambda;       // coefficients for function h(l,Z)
    double mu, Q0;            // coefficients for function Q(Z)
    double u1, u2, u3, u4;    // coefficients for function tau(Z)
    double cutsq;
    int ielement, jelement, kelement, dummy;    // dummy added for better alignment
  };

  double *preForceCoord;

  double cutmax;    // max cutoff for all elements
  Param *params;    // parameter set for an I-J-K interaction

  void allocate();
  void allocatePreLoops();
  void deallocatePreLoops();

  void read_file(char *);
  void setup() override;

  void edip_pair(double, double, const Param &, double &, double &, double &);
  void edip_fc(double, const Param &, double &, double &);
  void edip_fcut2(double, const Param &, double &, double &);
  void edip_tau(double, const Param &, double &, double &);
  void edip_h(double, double, const Param &, double &, double &, double &);
  void edip_fcut3(double, const Param &, double &, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif
