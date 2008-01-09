/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_WALL_GRAN_H
#define FIX_WALL_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallGran : public Fix {
 public:
  FixWallGran(class LAMMPS *, int, char **);
  ~FixWallGran();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void reset_dt();

 private:
  int wallstyle,pairstyle,wiggle,axis;
  double xkk,xkkt,gamman,xmu;
  double lo,hi,cylradius;
  double dt,gamman_dl,gammas_dl;
  double amplitude,period,omega,time_origin;

  int *touch;
  double **shear;

  void no_history(double, double, double, double, double *,
		  double *, double *, double *, double *, double, double);
  void history(double, double, double, double, double *,
	       double *, double *, double *, double *, double, double,
	       double *);
  void hertzian(double, double, double, double, double *,
		double *, double *, double *, double *, double, double,
		double *);
};

}

#endif
