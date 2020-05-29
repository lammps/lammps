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

FixStyle(rattle,FixRattle)

#else

#ifndef LMP_FIX_RATTLE_H
#define LMP_FIX_RATTLE_H

#include "fix_shake.h"

namespace LAMMPS_NS {

class FixRattle : public FixShake {
 public:
  double **vp;                // array for unconstrained velocities
  int comm_mode;              // mode for communication pack/unpack
  double derr_max;            // distance error
  double verr_max;            // velocity error

  FixRattle(class LAMMPS *, int, char **);
  ~FixRattle();
  int setmask();
  virtual void init();
  virtual void post_force(int);
  virtual void post_force_respa(int, int, int);
  virtual void final_integrate();
  virtual void final_integrate_respa(int,int);

  virtual void correct_coordinates(int vflag);
  virtual void correct_velocities();
  virtual void shake_end_of_step(int vflag);

  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  virtual void reset_dt();

 private:
  void update_v_half_nocons();
  void update_v_half_nocons_respa(int);

  void vrattle2(int m);
  void vrattle3(int m);
  void vrattle4(int m);
  void vrattle3angle(int m);
  void solve3x3exactly(const double a[][3], const double c[], double l[]);
  void solve2x2exactly(const double a[][2], const double c[], double l[]);

  // debugging methods

  bool check3angle(double ** v, int m, bool checkr, bool checkv);
  bool check2(double **v, int m, bool checkr, bool checkv);
  bool check3(double **v, int m, bool checkr, bool checkv);
  bool check4(double **v, int m, bool checkr, bool checkv);
  bool check_constraints(double **v, bool checkr, bool checkv);
  void end_of_step();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Fix rattle should come after all other integration fixes

UNDOCUMENTED

E: Rattle determinant = 0.0

The determinant of the matrix being solved for a single cluster
specified by the fix rattle command is numerically invalid.

E: Rattle failed

UNDOCUMENTED

E: Coordinate constraints are not satisfied up to desired tolerance

UNDOCUMENTED

E: Velocity constraints are not satisfied up to desired tolerance

UNDOCUMENTED

E: Velocity constraints are not satisfied up to desired tolerance!

UNDOCUMENTED

U: Fix rattle should come after all other integration fixes

This fix is designed to work after all other integration fixes change
atom positions.  Thus it should be the last integration fix specified.
If not, it will not satisfy the desired constraints as well as it
otherwise would.

U: Rattle failed

Certain constraints were not satisfied.

U: Coordinate constraints are not satisfied up to desired tolerance

Self-explanatory.

U: Rattle velocity constraints are not satisfied up to desired tolerance

Self-explanatory.

*/
