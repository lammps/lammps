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
FixStyle(rattle,FixRattle);
// clang-format on
#else

#ifndef LMP_FIX_RATTLE_H
#define LMP_FIX_RATTLE_H

#include "fix_shake.h"

namespace LAMMPS_NS {

class FixRattle : public FixShake {
 public:
  double **vp;        // array for unconstrained velocities
  int comm_mode;      // mode for communication pack/unpack
  double derr_max;    // distance error
  double verr_max;    // velocity error

  FixRattle(class LAMMPS *, int, char **);
  ~FixRattle() override;
  int setmask() override;
  void init() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void final_integrate() override;
  void final_integrate_respa(int, int) override;

  void correct_coordinates(int vflag) override;
  void correct_velocities() override;
  void shake_end_of_step(int vflag) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void reset_dt() override;

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

  bool check3angle(double **v, int m, bool checkr, bool checkv);
  bool check2(double **v, int m, bool checkr, bool checkv);
  bool check3(double **v, int m, bool checkr, bool checkv);
  bool check4(double **v, int m, bool checkr, bool checkv);
  bool check_constraints(double **v, bool checkr, bool checkv);
  void end_of_step() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
