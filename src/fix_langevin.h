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

FixStyle(langevin,FixLangevin)

#else

#ifndef LMP_FIX_LANGEVIN_H
#define LMP_FIX_LANGEVIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevin : public Fix {
 public:
  FixLangevin(class LAMMPS *, int, char **);
  virtual ~FixLangevin();
  int setmask();
  void init();
  void setup(int);
  virtual void initial_integrate(int);
  virtual void post_force(int);
  void post_force_respa(int, int, int);
  virtual void end_of_step();
  void reset_target(double);
  void reset_dt();
  int modify_param(int, char **);
  virtual double compute_scalar();
  double memory_usage();
  virtual void *extract(const char *, int &);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 protected:
  int gjfflag,nvalues,osflag,oflag,tallyflag,zeroflag,tbiasflag;
  int flangevin_allocated;
  double ascale;
  double t_start,t_stop,t_period,t_target;
  double *gfactor1,*gfactor2,*ratio;
  double energy,energy_onestep;
  double tsqrt;
  int tstyle,tvar;
  double gjfa, gjfsib; //gjf a and gjf sqrt inverse b
  char *tstr;

  class AtomVecEllipsoid *avec;

  int maxatom1,maxatom2;
  double **flangevin;
  double *tforce;
  double **franprev;
  double **lv; //half step velocity

  char *id_temp;
  class Compute *temperature;

  int nlevels_respa;
  class RanMars *random;
  int seed;

  template < int Tp_TSTYLEATOM, int Tp_GJF, int Tp_TALLY,
             int Tp_BIAS, int Tp_RMASS, int Tp_ZERO >
  void post_force_templated();

  void omega_thermostat();
  void angmom_thermostat();
  void compute_target();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix langevin period must be > 0.0

The time window for temperature relaxation must be > 0

E: Fix langevin omega requires atom style sphere

Self-explanatory.

E: Fix langevin angmom requires atom style ellipsoid

Self-explanatory.

E: Variable name for fix langevin does not exist

Self-explanatory.

E: Variable for fix langevin is invalid style

It must be an equal-style variable.

E: Fix langevin omega requires extended particles

One of the particles has radius 0.0.

E: Fix langevin angmom requires extended particles

This fix option cannot be used with point particles.

E: Cannot zero Langevin force of 0 atoms

The group has zero atoms, so you cannot request its force
be zeroed.

E: Fix langevin variable returned negative temperature

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

E: Fix langevin gjf cannot have period equal to dt/2

If the period is equal to dt/2 then division by zero will happen.

E: Fix langevin gjf should come before fix nve

Self-explanatory

E: Fix langevin gjf and respa are not compatible

Self-explanatory

W: Group for fix_modify temp != fix group

The fix_modify command is specifying a temperature computation that
computes a temperature on a different group of atoms than the fix
itself operates on.  This is probably not what you want to do.

*/
