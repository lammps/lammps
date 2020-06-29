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
  void final_integrate_respa(int, int);

  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);

  void reset_dt();

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
  int extra_flag,omega_flag,angmom_flag;
  int radius_flag,ellipsoid_flag,line_flag,tri_flag,body_flag;
  int theta_flag,quat_flag;
  int nlevels_respa,nrestart;
  int time_origin;

  double **xoriginal;         // original coords of atoms
  double *toriginal;          // original theta of atoms
  double **qoriginal;         // original quat of atoms
  int displaceflag,velocityflag;
  int maxatom;
  double **displace,**velocity;

  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix move cannot set linear z motion for 2d problem

Self-explanatory.

E: Fix move cannot set wiggle z motion for 2d problem

Self-explanatory.

E: Fix move cannot rotate around non z-axis for 2d problem

UNDOCUMENTED

E: Fix move cannot define z or vz variable for 2d problem

Self-explanatory.

E: Zero length rotation vector with fix move

Self-explanatory.

E: Variable name for fix move does not exist

Self-explanatory.

E: Variable for fix move is invalid style

Only equal-style variables can be used.

E: Cannot add atoms to fix move variable

Atoms can not be added afterwards to this fix option.

E: Resetting timestep size is not allowed with fix move

This is because fix move is moving atoms based on elapsed time.

U: Fix move cannot rotate aroung non z-axis for 2d problem

Self-explanatory.

*/
