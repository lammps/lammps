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

#ifdef INTEGRATE_CLASS

IntegrateStyle(respa,Respa)

#else

#ifndef LMP_RESPA_H
#define LMP_RESPA_H

#include "integrate.h"

namespace LAMMPS_NS {

class Respa : public Integrate {
 public:
                          // public so Fixes, Pairs, Neighbor can see them
  int nlevels;            // number of rRESPA levels
                          // 0 = innermost level, nlevels-1 = outermost level
  double *step;           // timestep at each level
  int *loop;              // sub-cycling factor at each level
  double cutoff[4];       // cutoff[0] and cutoff[1] = between inner and middle
                          // cutoff[2] and cutoff[3] = between middle and outer
                          // if no middle then 0,1 = 2,3

  int level_bond,level_angle,level_dihedral;   // level to compute forces at
  int level_improper,level_pair,level_kspace;
  int level_inner,level_middle,level_outer;

  Respa(class LAMMPS *, int, char **);
  virtual ~Respa();
  virtual void init();
  virtual void setup();
  virtual void setup_minimal(int);
  virtual void run(int);
  virtual void cleanup();
  virtual void reset_dt();

  void copy_f_flevel(int);
  void copy_flevel_f(int);

 protected:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int torqueflag,erforceflag;
  int e_flag,rho_flag;

  int *newton;                      // newton flag at each level
  class FixRespa *fix_respa;        // Fix to store the force level array

  virtual void recurse(int);
  void force_clear(int);
  void sum_flevel_f();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Respa levels must be >= 1

Self-explanatory.

E: Cannot set both respa pair and inner/middle/outer

In the rRESPA integrator, you must compute pairwise potentials either
all together (pair), or in pieces (inner/middle/outer).  You can't do
both.

E: Must set both respa inner and outer

Cannot use just the inner or outer option with respa without using the
other.

E: Cannot set respa middle without inner/outer

In the rRESPA integrator, you must define both a inner and outer
setting in order to use a middle setting.

E: Invalid order of forces within respa levels

For respa, ordering of force computations within respa levels must
obey certain rules.  E.g. bonds cannot be compute less frequently than
angles, pairwise forces cannot be computed less frequently than
kspace, etc.

W: One or more respa levels compute no forces

This is computationally inefficient.

E: Respa inner cutoffs are invalid

The first cutoff must be <= the second cutoff.

E: Respa middle cutoffs are invalid

The first cutoff must be <= the second cutoff.

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

E: Pair style does not support rRESPA inner/middle/outer

You are attempting to use rRESPA options with a pair style that
does not support them.

U: Fix shake with rRESPA computes invalid pressures

This is a known bug in LAMMPS that has not yet been fixed.  If you use
SHAKE with rRESPA and perform a constant volume simulation (e.g. using
fix npt) this only affects the output pressure, not the dynamics of
the simulation.  If you use SHAKE with rRESPA and perform a constant
pressure simulation (e.g. using fix npt) then you will be
equilibrating to the wrong volume.

*/
