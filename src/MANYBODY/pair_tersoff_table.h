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

/* ----------------------------------------------------------------------
   Contributing author: Luca Ferraro (CASPUR)
   email: luca.ferraro@caspur.it

   Tersoff Potential
   References:
    1) Tersoff, Phys. Rev. B 39, 5566 (1988)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tersoff/table,PairTersoffTable);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_TABLE_H
#define LMP_PAIR_TERSOFF_TABLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTersoffTable : public Pair {
 public:
  PairTersoffTable(class LAMMPS *);
  ~PairTersoffTable() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  static constexpr int NPARAMS_PER_LINE = 17;

 protected:
  struct Param {
    double lam1, lam2, lam3;
    double c, d, h;
    double gamma, powerm;
    double powern, beta;
    double biga, bigb, cutoffR, cutoffS;
    double cut, cutsq;
    int ielement, jelement, kelement;
    int powermint;
  };

  double cutmax;    // max cutoff for all elements
  Param *params;    // parameter set for an I-J-K interaction

  void allocate();

  void read_file(char *);
  void setup_params();

  // pre-loop coordination functions

  double **preGtetaFunction, **preGtetaFunctionDerived;
  double *preCutoffFunction, *preCutoffFunctionDerived;
  virtual void allocatePreLoops();
  virtual void deallocatePreLoops();

  // grids

  double minArgumentExponential;
  double *exponential, ***cutoffFunction, ***cutoffFunctionDerived;
  double **gtetaFunction, **gtetaFunctionDerived;
  double **betaZetaPower, **betaZetaPowerDerived;

  void allocateGrids();
  void deallocateGrids();
};

}    // namespace LAMMPS_NS

#endif
#endif
