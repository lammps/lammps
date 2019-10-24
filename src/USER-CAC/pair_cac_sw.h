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

PairStyle(cac/sw,PairCACSW)

#else

#ifndef LMP_PAIR_SW_CAC_H
#define LMP_PAIR_SW_CAC_H

#include "pair_cac.h"

namespace LAMMPS_NS {

class PairCACSW : public PairCAC {
 public:
  PairCACSW(class LAMMPS *);
  virtual ~PairCACSW();

  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);

	struct Param {
		double epsilon, sigma;
		double littlea, lambda, gamma, costheta;
		double biga, bigb;
		double powerp, powerq;
		double tol;
		double cut, cutsq;
		double sigma_gamma, lambda_epsilon, lambda_epsilon2;
		double c1, c2, c3, c4, c5, c6;
		int ielement, jelement, kelement;
	};

 protected:
	//potential params
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction

  double **cut;
 
  void allocate();
  void read_file(char *);
  virtual void setup_params();
  void twobody(Param *, double, double &, int, double &);
  void threebody(Param *, Param *, Param *, double, double, double *, double *,
  double *, double *, int, double &);
  //double density_map(double);
  
  void force_densities(int, double, double, double, double, double
    &fx, double &fy, double &fz);
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style cac/sw requires atom IDs

This is a requirement to use the SW potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Stillinger-Weber potential file %s

The specified SW potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Stillinger-Weber potential file

Incorrect number of words per line in the potential file.

E: Illegal Stillinger-Weber parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/