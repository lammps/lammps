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

FixStyle(poems,FixPOEMS)

#else

#ifndef LMP_FIX_POEMS_H
#define LMP_FIX_POEMS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPOEMS : public Fix  {
 public:
  FixPOEMS(class LAMMPS *, int narg, char **arg);
  ~FixPOEMS();
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void post_force(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void post_force_respa(int, int, int);
  void final_integrate_respa(int, int);

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double memory_usage();

  void pre_neighbor();
  int dof(int);
  void deform(int);
  int modify_param(int, char **);
  void reset_dt();

 private:
  int me;
  double dtv,dtf,dthalf;
  double *step_respa;
  int nlevels_respa;
  double total_ke;
  int earlyflag;    // 1 if forces and torques are computed at post_force()

  // atom assignment to rigid bodies
  // double count joint atoms as being in multiple bodies

  int *natom2body;         // # of bodies each atom is part of
  int **atom2body;         // list of bodies each atom is part of
  double **displace;       // atom displace in body coords for 1st body it's in

  // rigid body properties
  // only nrigid double counts joint atoms as being in multiple bodies
  // other quantities only count a joint atom as being in 1st body

  int nbody;                // # of rigid bodies
  int *nrigid;              // # of atoms in each rigid body
  double *masstotal;        // total mass of each rigid body
  double **xcm;             // coords of center-of-mass of each rigid body
  double **vcm;             // velocity of center-of-mass of each
  double **fcm;             // force on center-of-mass of each
  double **inertia;         // 6 inertia components of each (xx,yy,zz,xy,yz,xz)
  double **ex_space,**ey_space,**ez_space;
                            // orientation of each body in space coords
  double **angmom;          // angular momentum of each in space coords
  double **omega;           // angular velocity of each in space coords
  double **torque;          // torque on each rigid body in space coords
  double **sum,**all;       // work vectors

  // joint attributes between pairs of rigid bodies

  int ncluster;             // # of independent clusters of coupled bodies
  int njoint;               // # of interbody joints
  int **jointbody;          // indices of 2 rigid bodies in each joint (1-N)
  double **xjoint;          // coords of each joint point
  int nfree;                // # of isolated unconnected bodies
  int *freelist;            // indices of isolated bodies (1-N)

  // POEMS object

  class Workspace *poems;

  // internal class functions

  void compute_forces_and_torques();
  void readfile(char *);
  int readline(FILE *, char **, int *);
  void jointbuild();
  void sortlist(int, tagint **);
  int loopcheck(int, int, tagint **);
  int jacobi(double **, double *, double **);
  void rotate(double **, int, int, int, int, double, double);
  void omega_from_mq(double *, double *, double *, double *,
                     double *, double *);
  void set_v();
  void set_xv();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find fix poems group ID

A group ID used in the fix poems command does not exist.

E: Must use a molecular atom style with fix poems molecule

Self-explanatory.

E: Too many molecules for fix poems

The limit is 2^31 = ~2 billion molecules.

E: No rigid bodies defined

The fix specification did not end up defining any rigid bodies.

E: Atom in too many rigid bodies - boost MAXBODY

Fix poems has a parameter MAXBODY (in fix_poems.cpp) which determines
the maximum number of rigid bodies a single atom can belong to (i.e. a
multibody joint).  The bodies you have defined exceed this limit.

E: One or zero atoms in rigid body

Any rigid body defined by the fix rigid command must contain 2 or more
atoms.

W: More than one fix poems

It is not efficient to use fix poems more than once.

E: POEMS fix must come before NPT/NPH fix

NPT/NPH fix must be defined in input script after all poems fixes,
else the fix contribution to the pressure virial is incorrect.

E: Insufficient Jacobi rotations for POEMS body

Eigensolve for rigid body was not sufficiently accurate.

E: Rigid body has degenerate moment of inertia

Fix poems will only work with bodies (collections of atoms) that have
non-zero principal moments of inertia.  This means they must be 3 or
more non-collinear atoms, even with joint atoms removed.

E: Bad principal moments

Fix rigid did not compute the principal moments of inertia of a rigid
group of atoms correctly.

E: Cannot open fix poems file %s

The specified file cannot be opened.  Check that the path and name are
correct.

W: No joints between rigid bodies, use fix rigid instead

The bodies defined by fix poems are not connected by joints.  POEMS
will integrate the body motion, but it would be more efficient to use
fix rigid.

E: Cyclic loop in joint connections

Fix poems cannot (yet) work with coupled bodies whose joints connect
the bodies in a ring (or cycle).

E: Tree structure in joint connections

Fix poems cannot (yet) work with coupled bodies whose joints connect
the bodies in a tree structure.

*/
