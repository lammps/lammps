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
  virtual void init();
  void setup(int);
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);

  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  virtual void set_arrays(int);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(int, double *);
  virtual int pack_restart(int, double *);
  virtual void unpack_restart(int, int);
  virtual int size_restart(int);
  virtual int maxsize_restart();
  void reset_dt();

  void hooke(double, double, double, double, double *, double *,
             double *, double *, double *, double, double, double*);
  void hooke_history(double, double, double, double, double *,
                     double *, double *, double *, double *, double,
                     double, double *, double *);
  void hertz_history(double, double, double, double, double *,
                     double, double *, double *, double *, double *,
                     double, double, double *, double *);
  void granular(double, double, double, double, double *, double,
                double *, double *, double *, double *, double,
                double, double *, double *);

  double pulloff_distance(double);

 protected:
  int wallstyle,wiggle,wshear,axis;
  int pairstyle,nlevels_respa;
  bigint time_origin;
  double kn,kt,gamman,gammat,xmu;

  //For granular
  //Model choices
  int normal_model, damping_model;
  int tangential_model, roll_model, twist_model;

  int beyond_contact;

  //History flags
  int normal_history, tangential_history, roll_history, twist_history;

  //Indices of history entries
  int normal_history_index;
  int tangential_history_index;
  int roll_history_index;
  int twist_history_index;

  //Material coefficients
  double Emod, poiss, Gmod;

  //Contact model coefficients
  double normal_coeffs[4];
  double tangential_coeffs[3];
  double roll_coeffs[3];
  double twist_coeffs[3];

  double lo,hi,cylradius;
  double amplitude,period,omega,vshear;
  double dt;
  char *idregion;

  int use_history;       // if particle/wall interaction stores history
  int history_update;   // flag for whether shear history is updated
  int size_history;      // # of shear history values per contact

  // shear history for single contact per particle

  double **history_one;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // Store particle interactions
  int store;
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

E: Invalid fix wall/gran interaction style

UNDOCUMENTED

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/gran

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/gran

Self-explanatory.

E: Invalid shear direction for fix wall/gran

Self-explanatory.

E: Cannot wiggle or shear with fix wall/gran/region

UNDOCUMENTED

U: Fix wall/gran is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
