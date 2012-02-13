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

#ifdef FIX_CLASS

FixStyle(wall/gran,FixWallGran)

#else

#ifndef LMP_FIX_WALL_GRAN_H
#define LMP_FIX_WALL_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallGran : public Fix {
 public:
  FixWallGran(class LAMMPS *, int, char **);
  virtual ~FixWallGran();
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void reset_dt();

 protected:
  int wallstyle,pairstyle,wiggle,wshear,axis;
  double kn,kt,gamman,gammat,xmu;
  double lo,hi,cylradius;
  double amplitude,period,omega,vshear;
  double dt;
  int nlevels_respa;
  int time_origin;

  bigint laststep;
  int *touch;
  double **shear;
  int shearupdate;

  void hooke(double, double, double, double, double *,
	     double *, double *, double *, double *, double, double);
  void hooke_history(double, double, double, double, double *,
		     double *, double *, double *, double *, double, double,
		     double *);
  void hertz_history(double, double, double, double, double *,
		     double *, double *, double *, double *, double, double,
		     double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix wall/gran requires atom style sphere

Self-explanatory.

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/gran

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/gran

Self-explanatory.

E: Invalid shear direction for fix wall/gran

Self-explanatory.

E: Fix wall/gran is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
