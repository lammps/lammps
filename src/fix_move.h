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

FixStyle(move,FixMove)

#else

#ifndef LMP_FIX_MOVE_H
#define LMP_FIX_MOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMove : public Fix {
 public:
  FixMove(class LAMMPS *, int, char **);
  ~FixMove();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int);

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);

 private:
  char *xvarstr,*yvarstr,*zvarstr,*vxvarstr,*vyvarstr,*vzvarstr;
  int mstyle;
  int vxflag,vyflag,vzflag,axflag,ayflag,azflag;
  double vx,vy,vz,ax,ay,az;
  double period,omega_rotate;
  double point[3],axis[3],runit[3];
  double dt,dtv,dtf;
  int xvar,yvar,zvar,vxvar,vyvar,vzvar;
  int xvarstyle,yvarstyle,zvarstyle,vxvarstyle,vyvarstyle,vzvarstyle;
  int omega_flag,nlevels_respa;
  int time_origin;

  double **xoriginal;         // original coords of atoms
  int displaceflag,velocityflag;
  int maxatom;
  double **displace,**velocity;
};

}

#endif
#endif
