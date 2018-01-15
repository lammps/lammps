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

PairStyle(pair/spin/me,PairSpinMe)

#else

#ifndef LMP_PAIR_SPIN_ME_H
#define LMP_PAIR_SPIN_ME_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpinMe : public Pair {
 public:
  PairSpinMe(class LAMMPS *);
  virtual ~PairSpinMe();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_me(int, int, double [3], double [3], double [3], double [3]);  
  void compute_me_mech(int, int, double [3], double [3], double [3]);  
 
  int me_flag;				// me flag
  int me_mech_flag;			// mech calculation flag

  double cut_spin_me_global;		// global me cutoff
  double **cut_spin_me;			// me cutoff distance 

 protected:
  int newton_pair_spin; 
  double hbar;

  double **ME, **ME_mech;		// me coeff in eV
  double **v_mex, **v_mey, **v_mez;	// me direction

  void allocate();
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

E: Pair spin requires atom attributes sp, mumag

The atom style defined does not have these attributes.

*/
