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
PairStyle(exp6/rx,PairExp6rx);
// clang-format on
#else

#ifndef LMP_PAIR_EXP6_RX_H
#define LMP_PAIR_EXP6_RX_H

#include "pair.h"

namespace LAMMPS_NS {

class PairExp6rx : public Pair {
 public:
  PairExp6rx(class LAMMPS *);
  ~PairExp6rx() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  struct Param {
    double epsilon, rm, alpha;
    int ispecies;
    char *name, *potential;    // names of unique molecules and interaction type
    char *tablename;           // name of interaction table
    int potentialType;         // enumerated interaction potential type.
  };

 protected:
  enum { LINEAR };
  enum { NONE, EXPONENT, POLYNOMIAL };
  double cut_global;
  double **cut;
  double **epsilon, **rm, **alpha;
  double **rminv, **buck1, **buck2, **offset;

  virtual void allocate();
  int *mol2param;    // mapping from molecule to parameters
  int nparams;       // # of stored parameter sets
  int maxparam;      // max # of parameter sets
  Param *params;     // parameter set for an I-J-K interaction

  int nspecies;
  virtual void read_file(char *);
  void read_file2(char *);
  void setup() override;

  int isite1, isite2;
  char *site1, *site2;
  void getMixingWeights(int, double &, double &, double &, double &, double &, double &, double &,
                        double &, double &, double &, double &, double &, double &, double &,
                        double &, double &) const;
  double exponentR, exponentEpsilon;
  int scalingFlag;
  void exponentScaling(double, double &, double &) const;
  void polynomialScaling(double, double &, double &, double &) const;
  double *coeffAlpha, *coeffEps, *coeffRm;
  bool fractionalWeighting;

  inline double func_rin(const double &) const;
  inline double expValue(const double) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
