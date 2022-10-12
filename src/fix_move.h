/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(move,FixMove);
// clang-format on
#else

#ifndef LMP_FIX_MOVE_H
#define LMP_FIX_MOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMove : public Fix {
 public:
  FixMove(class LAMMPS *, int, char **);
  ~FixMove() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;

  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int maxsize_restart() override;
  int size_restart(int) override;

  void reset_dt() override;

 private:
  char *xvarstr, *yvarstr, *zvarstr, *vxvarstr, *vyvarstr, *vzvarstr;
  int mstyle;
  int vxflag, vyflag, vzflag, axflag, ayflag, azflag;
  double vx, vy, vz, ax, ay, az;
  double period, omega_rotate;
  double point[3], axis[3], runit[3];
  double dt, dtv, dtf;
  int xvar, yvar, zvar, vxvar, vyvar, vzvar;
  int xvarstyle, yvarstyle, zvarstyle, vxvarstyle, vyvarstyle, vzvarstyle;
  int extra_flag, omega_flag, angmom_flag;
  int radius_flag, ellipsoid_flag, line_flag, tri_flag, body_flag;
  int theta_flag, quat_flag;
  int nlevels_respa, nrestart;
  int time_origin;

  double **xoriginal;    // original coords of atoms
  double *toriginal;     // original theta of atoms
  double **qoriginal;    // original quat of atoms
  int displaceflag, velocityflag;
  int maxatom;
  double **displace, **velocity;

  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;
};

}    // namespace LAMMPS_NS

#endif
#endif
