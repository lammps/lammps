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

PairStyle(cac/buck,PairCACBuck)

#else

#ifndef LMP_PAIR_BUCK_CAC_H
#define LMP_PAIR_BUCK_CAC_H

#include "pair_cac.h"

namespace LAMMPS_NS {

class PairCACBuck : public PairCAC {
 public:
  PairCACBuck(class LAMMPS *);
  virtual ~PairCACBuck();

  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);

 protected:

  double **cut;
  double **a, **rho, **c;
  double **rhoinv, **buck1, **buck2, **offset;
  
  void allocate();
  void force_densities(int, double, double, double, double, double
  &fx, double &fy, double &fz);
  virtual void settings(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Unexpected argument in pair cac/buck invocation; only accepts cutoff and the 'one' keyword

Self-explanatory.  Check the input script. See the documentation for the proper syntax.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
