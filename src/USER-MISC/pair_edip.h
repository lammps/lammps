/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(edip,PairEDIP)

#else

#ifndef LMP_PAIR_EDIP_H
#define LMP_PAIR_EDIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEDIP : public Pair {
 public:
  PairEDIP(class LAMMPS *);
  virtual ~PairEDIP();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 protected:
  struct Param {
    double A, B;
    double cutoffA, cutoffC, cutsq;
    double alpha, beta;
    double eta, gamm, lambda, mu, rho, sigma, Q0;
    double u1, u2, u3, u4;
    int ielement,jelement,kelement;
  };

  double *preInvR_ij;
  double *preExp3B_ij;
  double *preExp3BDerived_ij;
  double *preExp2B_ij;
  double *preExp2BDerived_ij;
  double *prePow2B_ij;
  double *preForceCoord;

  // grids

  double *cutoffFunction;
  double *cutoffFunctionDerived;
  double *pow2B;
  double *exp2B;
  double *exp3B;
  double *qFunctionGrid;
  double *expMinusBetaZeta_iZeta_iGrid;
  double *tauFunctionGrid;
  double *tauFunctionDerivedGrid;

  // this should be removed for multi species parametrizations
  // since these parameters should be addressed through indexes
  // see also the PairEDIP::setup()

  double A;
  double B;
  double rho;
  double cutoffA;
  double cutoffC;
  double sigma;
  double lambda;
  double gamm;
  double eta;
  double Q0;
  double mu;
  double beta;
  double alpha;
  double u1;
  double u2;
  double u3;
  double u4;

  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  void allocate();
  void allocatePreLoops(void);
  void deallocatePreLoops(void);
  void allocateGrids(void);
  void deallocateGrids(void);
  void initGrids(void);

  void read_file(char *);
  void setup();
};

}

#endif
#endif
