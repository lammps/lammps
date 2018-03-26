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

PairStyle(pair/spin/soc/neel,PairSpinSocNeel)

#else

#ifndef LMP_PAIR_SPIN_SOC_NEEL_H
#define LMP_PAIR_SPIN_SOC_NEEL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpinSocNeel : public Pair {
 public:
  PairSpinSocNeel(class LAMMPS *);
  virtual ~PairSpinSocNeel();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_soc_neel(int, int, double, double [3], double [3], double [3], double [3]);
  void compute_soc_mech_neel(int, int, double, double [3], double [3], double [3], double [3]);
 
  int soc_neel_flag;		// soc neel flag
  int soc_mech_flag; 		// mech calculation flag

  double cut_soc_global;
  double **cut_soc_neel;	// cutoff distance exchange

 protected:
  double hbar;

  // pseudo-dipolar coeff.
  double **g1, **g1_mech; 	// exchange coeffs gij
  double **g2, **g3; 		// g1 in eV, g2 adim, g3 in Ang

  // pseudo-quadrupolar coeff.
  double **q1, **q1_mech; 	// exchange coeffs qij
  double **q2, **q3; 		// q1 in eV, q2 adim, q3 in Ang

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

E: Pair spin requires atom attribute spin

The atom style defined does not have these attributes.

*/
