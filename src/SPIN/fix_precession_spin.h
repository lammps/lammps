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

FixStyle(precession/spin,FixPrecessionSpin)

#else

#ifndef LMP_FIX_PRECESSION_SPIN_H
#define LMP_FIX_PRECESSION_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPrecessionSpin : public Fix {
  friend class FixPour;

 public:
  FixPrecessionSpin(class LAMMPS *, int, char **);
  ~FixPrecessionSpin();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

  int zeeman_flag, aniso_flag, cubic_flag;
  void compute_single_precession(int, double *, double *);
  void compute_zeeman(int, double *);
  
  // uniaxial aniso calculations

  void compute_anisotropy(double *, double *);
  double compute_anisotropy_energy(double *);

  // cubic aniso calculations

  void compute_cubic(double *, double *);
  double compute_cubic_energy(double *);
 
 protected:
  int style; 			// style of the magnetic precession

  double degree2rad;
  double hbar;
  int ilevel_respa;
  int time_origin;
  int eflag;
  double eprec, eprec_all;

  int varflag;
  int magfieldstyle;
  int magvar;
  char *magstr;

  // zeeman field intensity and direction

  double H_field;
  double nhx, nhy, nhz;
  double hx, hy, hz; 		// temp. force variables

  // magnetic anisotropy intensity and direction

  double Ka;			// aniso const. in eV
  double Kah;			// aniso const. in rad.THz
  double nax, nay, naz;
  double Kax, Kay, Kaz; 	// temp. force variables

  // cubic anisotropy intensity

  double k1c,k2c;		// cubic const. in eV
  double k1ch,k2ch;		// cubic const. in rad.THz
  double nc1x,nc1y,nc1z;
  double nc2x,nc2y,nc2z;
  double nc3x,nc3y,nc3z;

  void set_magneticprecession();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal precession/spin command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

precession/spin fix command has 7 arguments:
fix  ID  group  precession/spin  magnitude (T or eV)  style (zeeman or anisotropy)
direction (3 cartesian coordinates)
*/
