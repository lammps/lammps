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

PairStyle(exp6/rx,PairExp6rx)

#else

#ifndef LMP_PAIR_EXP6_RX_H
#define LMP_PAIR_EXP6_RX_H

#include "pair.h"

namespace LAMMPS_NS {

class PairExp6rx : public Pair {
 public:
  PairExp6rx(class LAMMPS *);
  virtual ~PairExp6rx();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  struct Param {
    double epsilon,rm,alpha;
    int ispecies;
    char *name, *potential;      // names of unique molecules and interaction type
    char *tablename;             // name of interaction table
   int potentialType;              // enumerated interaction potential type.
  };

 protected:
  enum{LINEAR};
  enum{NONE,EXPONENT,POLYNOMIAL};
  double cut_global;
  double **cut;
  double **epsilon,**rm,**alpha;
  double **rminv,**buck1,**buck2,**offset;

  virtual void allocate();
  int *mol2param;               // mapping from molecule to parameters
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  int nspecies;
  virtual void read_file(char *);
  void read_file2(char *);
  void setup();

  int isite1, isite2;
  char *site1, *site2;
  void getMixingWeights(int, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &, double &) const;
  double exponentR, exponentEpsilon;
  int scalingFlag;
  void exponentScaling(double, double &, double &) const;
  void polynomialScaling(double, double &, double &, double &) const;
  double *coeffAlpha, *coeffEps, *coeffRm;
  bool fractionalWeighting;

  inline double func_rin(const double &) const;
  inline double expValue(const double) const;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E:  alpha_ij is 6.0 in pair exp6

Self-explanatory

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: PairExp6rx requires a fix rx command

The fix rx command must come before the pair style command in the input file

E:  There are no rx species specified

There must be at least one species specified through the fix rx command

E:  Site1 name not recognized in pair coefficients

The site1 keyword does not match the species keywords specified throug the fix rx command

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E:  Cannot open exp6/rx potential file %s

Self-explanatory

E:  Incorrect format in exp6/rx potential file

Self-explanatory

E:  Illegal exp6/rx parameters.  Rm and Epsilon must be greater than zero.  Alpha cannot be negative.

Self-explanatory

E:  Illegal exp6/rx parameters.  Interaction potential does not exist.

Self-explanatory

E:  Potential file has duplicate entry.

Self-explanatory

E:  The number of molecules in CG particle is less than 10*DBL_EPSILON.

Self-explanatory.  Check the species concentrations have been properly set
and check the reaction kinetic solver parameters in fix rx to more for
sufficient accuracy.


*/
