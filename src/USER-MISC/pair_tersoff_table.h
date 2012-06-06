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

/* ----------------------------------------------------------------------
   Contributing author: Luca Ferraro (CASPUR)
   email: luca.ferraro@caspur.it

   Tersoff Potential
   References:
    1) Tersoff, Phys. Rev. B 39, 5566 (1988)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tersoff/table,PairTersoffTable)

#else

#ifndef LMP_PAIR_Tersoff_H
#define LMP_PAIR_Tersoff_H

#include "pair.h"

namespace LAMMPS_NS {

class PairTersoffTable : public Pair {
 public:
  PairTersoffTable(class LAMMPS *);
  virtual ~PairTersoffTable();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 protected:
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,cutoffR, cutoffS;
    double cut,cutsq;
    int ielement,jelement,kelement;
    int powermint;
  };


  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  void allocate();

  void read_file(char *);
  void setup();

  // pre-loop coordination functions

  double **preGtetaFunction, **preGtetaFunctionDerived;
  double *preCutoffFunction, *preCutoffFunctionDerived;
  void allocatePreLoops(void);
  void deallocatePreLoops(void);

  // grids

  double minArgumentExponential;
  double *exponential, ***cutoffFunction, ***cutoffFunctionDerived;
  double **gtetaFunction, **gtetaFunctionDerived;
  double **betaZetaPower, **betaZetaPowerDerived;

  void allocateGrids(void);
  void deallocateGrids(void);
};

}

#endif
#endif
