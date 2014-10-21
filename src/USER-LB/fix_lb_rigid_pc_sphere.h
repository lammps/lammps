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

FixStyle(lb/rigid/pc/sphere,FixLbRigidPCSphere)

#else

#ifndef LMP_FIX_LB_RIGID_PC_SPHERE_H
#define LMP_FIX_LB_RIGID_PC_SPHERE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLbRigidPCSphere : public Fix {
 public:
  FixLbRigidPCSphere(class LAMMPS *, int, char **);
  virtual ~FixLbRigidPCSphere();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual double compute_scalar();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void pre_neighbor();
  int dof(int);
  void reset_dt();
  double compute_array(int, int);


 private:
  double **up;
  double **up_old;
  double *Gamma_MD;
  double expminusdttimesgamma,DMDcoeff;
  double expminusdttimesgammadiv2;
  double force_factor,torque_factor;

  double dtv,dtf;

  int nbody;                // # of rigid bodies
  int *nrigid;              // # of atoms in each rigid body
  int *nrigid_shell;
  double *masstotal;        // total mass of each rigid body
  double *masstotal_shell;
  double *sphereradius;
  double **xcm;             // coords of center-of-mass of each rigid body
  double **xcm_old;
  double **vcm;             // velocity of center-of-mass of each
  double **ucm;
  double **ucm_old;
  double **fcm;             // force on center-of-mass of each
  double **fcm_old;
  double **fcm_fluid;
  double **omega;          // angular momentum of each in space coords
  double **torque;          // torque on each rigid body in space coords
  double **torque_old;
  double **torque_fluid;
  double **torque_fluid_old;
  double **rotate;
  imageint *imagebody;               // image flags of xcm of each rigid body
  double **fflag;           // flag for on/off of center-of-mass force
  double **tflag;           // flag for on/off of center-of-mass torque

  int *body;                // which body each atom is part of (-1 if none)

  double **sum,**all;       // work vectors for each rigid body
  int **remapflag;          // PBC remap flags for each rigid body

  double tfactor;           // scale factor on temperature of rigid bodies

  int inner_nodes;          // ==1 if certain particle are inside the rigid
                            //  body and should not interact with the fluid. 
                            //  ==0 otherwise.
  int igroupinner;          // specifies the particles which are inside the
                            //  spherical rigid body, and do not interact with
                            //  the fluid.

  void set_xv();
  void set_v();
  
  void compute_up();

  class FixLbFluid *fix_lb_fluid;
};

}

#endif
#endif
