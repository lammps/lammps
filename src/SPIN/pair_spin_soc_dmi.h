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

PairStyle(pair/spin/soc/dmi,PairSpinSocDmi)

#else

#ifndef LMP_PAIR_SPIN_SOC_DMI_H
#define LMP_PAIR_SPIN_SOC_DMI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairSpinSocDmi : public Pair {
 public:
  PairSpinSocDmi(class LAMMPS *);
  virtual ~PairSpinSocDmi();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  void compute_dmi(int, int, double *, double *, double *, double *);
 
  int dmi_flag;    // dmi flag

  double cut_spin_dmi_global;    // short range pair cutoff
  double **cut_spin_dmi;      // cutoff distance dmi

 protected:
  int newton_pair_spin; 
  double hbar;

  double **DM;                     // dmi coeff in eV 
  double **v_dmx, **v_dmy, **v_dmz;// dmi direction

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
