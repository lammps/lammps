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
FixStyle(wall/gran,FixWallGran);
// clang-format on
#else

#ifndef LMP_FIX_WALL_GRAN_H
#define LMP_FIX_WALL_GRAN_H

#include "granular_model.h"
#include "fix.h"

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class FixWallGran : public Fix {
 public:
  FixWallGran(class LAMMPS *, int, char **);
  ~FixWallGran() override;
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

 protected:
  int wallstyle, wiggle, wshear, axis;
  int nlevels_respa;
  bigint time_origin;

  // for granular model choices
  class Granular_NS::GranularModel *model;

  double lo, hi, cylradius;
  double amplitude, period, omega, vshear;
  double dt;
  double Twall;
  char *idregion;

  int use_history;       // if particle/wall interaction stores history
  int history_update;    // flag for whether shear history is updated
  int size_history;      // # of shear history values per contact
  int heat_flag;

  int tvar;
  char *tstr;

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
