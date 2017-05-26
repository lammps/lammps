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

PairStyle(pair/spin,PairSpin)

#else

#ifndef LMP_PAIR_SPIN_H
#define LMP_PAIR_SPIN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpin : public Pair {
 public:
  PairSpin(class LAMMPS *);
  virtual ~PairSpin();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  //Test transf. force
//#define TRANS_FORCE
#if defined TRANS_FORCE
  void transferfm(double **);
#endif
  
 protected:
  double cut_spin_exchange_global, cut_spin_dipolar_global; //Global cutting distance
  double **cut_spin_exchange; //cutting distance for each exchange interaction
  double **cut_spin_dipolar;  //cutting distance for the dipolar interaction
  
  double **J_1, **J_2, **J_3; //coefficients for computing the exchange interaction Jij
                              //J1 is an energy (in eV), J2 is adim and J3 is a distance (in Ang)

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_style command

Self-explanatory.

E: Cannot (yet) use 'electron' units with spins

This feature is not yet supported.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair spin requires atom attributes sp, mumag

The atom style defined does not have these attributes.

*/
