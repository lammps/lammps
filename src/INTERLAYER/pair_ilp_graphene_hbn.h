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
PairStyle(ilp/graphene/hbn,PairILPGrapheneHBN);
// clang-format on
#else

#ifndef LMP_PAIR_ILP_GRAPHENE_HBN_H
#define LMP_PAIR_ILP_GRAPHENE_HBN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBN : public Pair {
 public:
  PairILPGrapheneHBN(class LAMMPS *);
  ~PairILPGrapheneHBN() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void calc_FvdW(int, int);
  double single(int, int, int, int, double, double, double, double &) override;

  static constexpr int NPARAMS_PER_LINE = 13;

  enum { ILP_GrhBN, ILP_TMD, SAIP_METAL };    // for telling class variants apart in shared code

 protected:
  int me;
  int variant;
  int maxlocal;            // size of numneigh, firstneigh arrays
  int pgsize;              // size of neighbor page
  int oneatom;             // max # of neighbors for one atom
  MyPage<int> *ipage;      // neighbor list pages
  int *ILP_numneigh;       // # of pair neighbors for each atom
  int **ILP_firstneigh;    // ptr to 1st neighbor of each atom
  int tap_flag;            // flag to turn on/off taper function

  struct Param {
    double z0, alpha, epsilon, C, delta, d, sR, reff, C6, S;
    double delta2inv, seff, lambda, rcut;
    int ielement, jelement;
  };
  Param *params;    // parameter set for I-J interactions
  int nmax;         // max # of atoms

  double cut_global;
  double cut_normal;
  double **cutILPsq;    // mapping the cutoff from element pairs to parameters
  double **offset;
  double **normal;
  double ***dnormdri;
  double ****dnormal;

  // adds for ilp/tmd
  int Nnei;    // max # of nearest neighbors for one atom
  double **dnn;
  double **vect;
  double **pvet;
  double ***dpvet1;
  double ***dpvet2;
  double ***dNave;

  virtual void ILP_neigh();
  virtual void calc_normal();
  virtual void calc_FRep(int, int);
  void read_file(char *);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
