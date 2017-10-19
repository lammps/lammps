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

PairStyle(pair/spin/soc/landau,PairSpinSocLandau)

#else

#ifndef LMP_PAIR_SPIN_SOC_LANDAU_H
#define LMP_PAIR_SPIN_SOC_LANDAU_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpinSocLandau : public Pair {
 public:
  PairSpinSocLandau(class LAMMPS *);
  virtual ~PairSpinSocLandau();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_soc_neel(int, int, double, double *, double *, double *,double *, double *);
  void compute_soc_mech_neel(int, int, double, double *, double *, double *,double *, double *);
 
  int soc_neel_flag; // soc neel flag
  int mech_flag; // mech calc. flag

  double cut_soc_global;
  double **cut_soc_neel; // cutoff distance exchange

 protected:
  int newton_pair_spin; 
  double hbar;

  double **K1, **K1_mech; // exchange coeffs Kij
  double **K2, **K3; // K1 in eV, K2 adim, K3 in Ang

  double *spi, *spj; // temp. spin vals. in compute
  double *fi, *fj;   // temp. mech. forces  in compute
  double *fmi, *fmj; // temp. mag. forces in compute
  double *rij;       // norm. inter atomic vectors

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
