/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(WALL/GRAN/OLD,FixWallGranOld);
// clang-format on
#else

#ifndef LMP_FIX_WALL_GRAN_OLD_H
#define LMP_FIX_WALL_GRAN_OLD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallGranOld : public Fix {
 public:
  enum { HOOKE, HOOKE_HISTORY, HERTZ_HISTORY, GRANULAR };
  enum { NORMAL_NONE, NORMAL_HOOKE, NORMAL_HERTZ, HERTZ_MATERIAL, DMT, JKR };

  FixWallGranOld(class LAMMPS *, int, char **);
  ~FixWallGranOld() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  void reset_dt() override;

  void hooke(double, double, double, double, double *, double *, double *, double *, double *,
             double, double, double *);
  void hooke_history(double, double, double, double, double *, double *, double *, double *,
                     double *, double, double, double *, double *);
  void hertz_history(double, double, double, double, double *, double, double *, double *, double *,
                     double *, double, double, double *, double *);
  void granular(double, double, double, double, double *, double, double *, double *, double *,
                double *, double, double, double *, double *);

  double pulloff_distance(double);

 protected:
  int wallstyle, wiggle, wshear, axis;
  int pairstyle, nlevels_respa;
  bigint time_origin;
  double kn, kt, gamman, gammat, xmu;

  // for granular model choices
  int normal_model, damping_model;
  int tangential_model, roll_model, twist_model;
  int limit_damping;

  // history flags
  int normal_history, tangential_history, roll_history, twist_history;

  // indices of history entries
  int normal_history_index;
  int tangential_history_index;
  int roll_history_index;
  int twist_history_index;

  // material coefficients
  double Emod, poiss, Gmod;

  // contact model coefficients
  double normal_coeffs[4];
  double tangential_coeffs[3];
  double roll_coeffs[3];
  double twist_coeffs[3];

  double lo, hi, cylradius;
  double amplitude, period, omega, vshear;
  double dt;
  char *idregion;

  int use_history;       // if particle/wall interaction stores history
  int history_update;    // flag for whether shear history is updated
  int size_history;      // # of shear history values per contact

  // shear history for single contact per particle

  double **history_one;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, null pointer if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // store particle interactions

  int store;

  void clear_stored_contacts();
};

}    // namespace LAMMPS_NS

#endif
#endif
