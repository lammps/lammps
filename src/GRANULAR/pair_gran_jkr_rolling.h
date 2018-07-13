/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(gran/jkr/rolling,PairGranJKRRolling)

#else

#ifndef LMP_PAIR_GRAN_JKR_ROLLING_H
#define LMP_PAIR_GRAN_JKR_ROLLING_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranJKRRolling : public PairGranHookeHistory {
public:
  PairGranJKRRolling(class LAMMPS *);
  virtual ~PairGranJKRRolling();
  virtual void compute(int, int);
  void settings(int, char **); //Eventually set this through coeff method so that user can specify a particular i-j set of coefficients
  double single(int, int, int, int, double, double, double, double &);
  double *E_one, *G_one, *pois, *muS_one, *cor, *alpha_one, *Ecoh_one, *kR_one, *muR_one, *etaR_one; //Public so as to be accessible to fix/wall/gran
private:
  double **E, **G, **alpha, **muS, **Ecoh, **kR, **muR, **etaR, **gamman;
  int normaldamp, rollingdamp;



};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

 */
