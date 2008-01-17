/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#ifndef FIX_POEMS_H
#define FIX_POEMS_H

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
  void final_integrate_respa(int);

  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double memory_usage();

  void pre_neighbor();
  int dof(int);
  void deform(int);
  void reset_dt();

 private:
  int me;
  double dtv,dtf,dthalf;
  double *step_respa;
  int nlevels_respa;
  double total_ke;

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

  void readfile(char *);
  int readline(FILE *, char **, int *);
  void jointbuild();
  void sortlist(int, int **);
  int loopcheck(int, int, int **);
  int jacobi(double **, double *, double **);
  void rotate(double **, int, int, int, int, double, double);
  void omega_from_mq(double *, double *, double *, double *,
		     double *, double *);
  void set_v();
  void set_xv();
};

}

#endif
