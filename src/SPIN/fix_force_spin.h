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

 protected:
  int style;
  
  double xmag, ymag, zmag; //Magnetic force
  double degree2rad;
  int ilevel_respa;
  int time_origin;
  int eflag;
  double emag, emag_all;

  int varflag;
  int magfieldstyle;
  int magvar;
  char *magstr;
   
  double H_field; //Zeeman field intensity and direction
  double Hx, Hy, Hz;
  
  double Ka; //Magnetic anisotropy intensity and direction
  double Kax, Kay, Kaz;

  void set_magneticforce();
   
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix force/spin does not exist

Self-explanatory.

E: Variable for fix force/spin is invalid style

Only equal-style variables can be used.

*/
