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
  virtual void compute_magnetomech(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_exchange(int, int, double, double *, double *,double *, double *);
  void compute_exchange_mech(int, int, double, double *, double *, double *,double *, double *);
  void compute_dmi(int, int, double *, double *, double *, double *);
  void compute_me(int, int, double *, double *, double *, double *);  
 
 //Test for seq. integ.
 //protected: 
  int exch_flag,dmi_flag,me_flag;
  double cut_spin_pair_global;
  double cut_spin_dipolar_global;
  
  double **cut_spin_exchange; // cutoff distance exchange
  double **cut_spin_dmi;      // cutoff distance dmi
  double **cut_spin_me;       // cutoff distance me 

 protected: 
  double **J1_mag, **J1_mech; // exchange coeffs Jij
  double **J2, **J3; // J1 in eV, J2 adim, J3 in Ang

  double **DM;                     // dmi coeff in eV 
  double **v_dmx, **v_dmy, **v_dmz;// dmi direction

  double **ME;                     // me coeff in eV
  double **v_mex, **v_mey, **v_mez;// me direction

  double *spi, *spj; // temp. spin vals. in compute
  double *fi, *fj;   // temp. mech. forces  in compute
  double *fmi, *fmj; // temp. mag. forces in compute
  double *del;       // inter atomic distances

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
