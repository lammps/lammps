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
FixStyle(settorque/atom,FixSetTorqueAtom);
// clang-format on
#else

#ifndef LMP_FIX_SET_TORQUE_ATOM_H
#define LMP_FIX_SET_TORQUE_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSetTorqueAtom : public Fix {
 public:
  FixSetTorqueAtom(class LAMMPS *, int, char **);
  ~FixSetTorqueAtom() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_vector(int) override;

  double memory_usage() override;

 protected:
  double xvalue, yvalue, zvalue;
  int varflag;
  char *xstr, *ystr, *zstr;
  char *idregion;
  class Region *region;
  int xvar, yvar, zvar, xstyle, ystyle, zstyle;
  double toriginal[3], toriginal_all[3], toriginal_saved[3];
  int torque_flag;
  int nlevels_respa, ilevel_respa;

  int maxatom;
  double **storque;
};

}    // namespace LAMMPS_NS

#endif
#endif
