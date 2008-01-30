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

#ifndef PAIR_SW_H
#define PAIR_SW_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSW : public Pair {
 public:
  PairSW(class LAMMPS *);
  ~PairSW();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 private:
  struct Param {
    double epsilon,sigma;
    double littlea,lambda,gamma,costheta;
    double biga,bigb;
    double powerp,powerq;
    double tol;
    double cut,cutsq;
    double sigma_gamma,lambda_epsilon,lambda_epsilon2;
    double c1,c2,c3,c4,c5,c6;
    int ielement,jelement,kelement;
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
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double *, double *,
		 double *, double *, int, double &);
};

}

#endif
