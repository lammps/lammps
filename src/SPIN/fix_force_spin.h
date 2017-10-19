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

#ifdef FIX_CLASS

FixStyle(force/spin,FixForceSpin)

#else

#ifndef LMP_FIX_FORCE_SPIN_H
#define LMP_FIX_FORCE_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixForceSpin : public Fix {
  friend class FixPour;

 public:
  FixForceSpin(class LAMMPS *, int, char **);
  ~FixForceSpin();
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  double compute_scalar();

  int zeeman_flag, aniso_flag;	
  void compute_zeeman(int, double *);
  void compute_anisotropy(int, double *, double *);

 protected:
  int style; // style of the magnetic force
  
  double degree2rad;
  double hbar;
  int ilevel_respa;
  int time_origin;
  int eflag;
  double emag;

  int varflag;
  int magfieldstyle;
  int magvar;
  char *magstr;
   
  // zeeman field intensity and direction
  double H_field; 
  double nhx, nhy, nhz;
  double hx, hy, hz; // temp. force variables
  
  // magnetic anisotropy intensity and direction
  double Ka; 
  double nax, nay, naz;
  double Kax, Kay, Kaz; // temp. force variables

  // temp. spin variables
  double *spi, *fmi;
  
  void set_magneticforce();
 
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal force/spin command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
