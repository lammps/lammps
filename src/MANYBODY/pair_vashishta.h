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
PairStyle(vashishta,PairVashishta);
// clang-format on
#else

#ifndef LMP_PAIR_VASHISHITA_H
#define LMP_PAIR_VASHISHITA_H

#include "pair.h"

namespace LAMMPS_NS {

class PairVashishta : public Pair {
 public:
  PairVashishta(class LAMMPS *);
  ~PairVashishta() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;

  static constexpr int NPARAMS_PER_LINE = 17;

  struct Param {
    double bigb, gamma, r0, bigc, costheta;
    double bigh, eta, zi, zj;
    double lambda1, bigd, mbigd, lambda4, bigw, cut;
    double lam1inv, lam4inv, zizj, heta, big2b, big6w;
    double rcinv, rc2inv, rc4inv, rc6inv, rceta;
    double cutsq2, cutsq;
    double lam1rc, lam4rc, vrcc2, vrcc3, vrc, dvrc, c0;
    int ielement, jelement, kelement;
  };

 protected:
  double cutmax;      // max cutoff for all elements
  Param *params;      // parameter set for an I-J-K interaction
  double r0max;       // largest value of r0
  int maxshort;       // size of short neighbor list array
  int *neighshort;    // short neighbor list array

  void allocate();
  void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double *, double *, double *, double *,
                 int, double &);
};

}    // namespace LAMMPS_NS

#endif
#endif
