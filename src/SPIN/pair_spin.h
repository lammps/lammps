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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)
   
   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018). 
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics 
   and molecular dynamics. arXiv preprint arXiv:1801.10233.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(pair/spin,PairSpin)

#else

#ifndef LMP_PAIR_SPIN_H
#define LMP_PAIR_SPIN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpin : public Pair {
friend class FixNVESpin;
 public:
  PairSpin(class LAMMPS *);
  virtual ~PairSpin();
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int) {return 0.0;}
  virtual void compute(int, int);


  // functions called from fix/nve/spin
  // used to evaluate force on spin i inside integration loops

  int init_pair();					// init call of PairSpin styles
  double larger_cutoff;					// larger Pair/Spin style cutoff
  double larger_spin_cutoff();				// compute larger_cutoff
  void compute_pair_single_spin(int, double *);		// compute pairs for one atom

  class Pair *pair;					// unused try class
  int nspinstyle;					// # of magnetic pairs
  class Pair *pair_spin_match(const char *, int, int);	// initializing nspinstyle
  char *pair_keyword;					// main pair style
  char **pair_spin_keywords;				// list of Pair/Spin style names

  // # and lists of Pair/Spin style classes

  int nexchangespinstyle;				// # of exchange pairs
  class PairSpinExchange **exchange_spin_styles; 	// list of Pair/Spin/Exchange2 classes

  int ndmispinstyle;					// # of dmi pairs
  class PairSpinDmi **dmi_spin_styles; 			// list of Pair/Spin/Dmi2 classes

  int nneelspinstyle;					// # of dmi pairs
  class PairSpinNeel **neel_spin_styles; 		// list of Pair/Spin/Dmi2 classes

  int nmespinstyle;					// # of me pairs
  class PairSpinMe **me_spin_styles; 			// list of Pair/Spin/Me2 classes

 protected:
  int lattice_flag; 				// if 0 spins only, 1 if spin-lattice
  double hbar;					// Planck constant (eV.ps.rad-1) 
  class FixNVESpin *lockfixnvespin;		// ptr to FixNVESpin2 for setups

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_spin command

Self-explanatory.

E: Spin simulations require metal unit style

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair spin requires atom attribute spin

The atom style defined does not have these attributes.

*/
