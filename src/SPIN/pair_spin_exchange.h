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

PairStyle(pair/spin/exchange,PairSpinExchange)

#else

#ifndef LMP_PAIR_SPIN_EXCHANGE_H
#define LMP_PAIR_SPIN_EXCHANGE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpinExchange : public Pair {
 public:
  PairSpinExchange(class LAMMPS *);
  virtual ~PairSpinExchange();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_exchange(int, int, double, double fmi[3], double spi[3], double spj[3]);
  void compute_exchange_mech(int, int, double, double rij[3], double fi[3], double spi[3], double spj[3]);
 
  int exch_flag;			// magnetic exchange flag
  int exch_mech_flag;			// mechanic exchange flags
  double cut_spin_exchange_global;	// global exchange cutoff
  double **cut_spin_exchange;		// cutoff distance per exchange

 protected:
  double hbar;
  double **J1_mag, **J1_mech;		// exchange coeffs Jij
  double **J2, **J3;			// J1 in eV, J2 adim, J3 in Ang

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
