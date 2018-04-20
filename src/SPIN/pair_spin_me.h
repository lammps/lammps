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

PairStyle(spin/me,PairSpinMe)

#else

#ifndef LMP_PAIR_SPIN_ME_H
#define LMP_PAIR_SPIN_ME_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinMe : public PairSpin {
 public:
  PairSpinMe(class LAMMPS *);
  virtual ~PairSpinMe();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_me(int, int, double, double *, double *, double *, double *);  
  void compute_me_mech(int, int, double *, double *, double *);  
 
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  
  double cut_spin_me_global;		// global me cutoff

 protected:
  double **ME, **ME_mech;		// me coeff in eV
  double **v_mex, **v_mey, **v_mez;	// me direction
  double **cut_spin_me;			// me cutoff distance 

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
