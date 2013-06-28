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

PairStyle(nb3b/harmonic,PairNb3bHarmonic)

#else

#ifndef LMP_PAIR_NB3B_HARMONIC_H
#define LMP_PAIR_NB3B_HARMONIC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairNb3bHarmonic : public Pair {
 public:
  PairNb3bHarmonic(class LAMMPS *);
  ~PairNb3bHarmonic();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

 private:
  struct Param {
    double k_theta, theta0, cutoff;
    double cut,cutsq;
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
#endif
