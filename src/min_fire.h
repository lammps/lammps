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

#ifdef MINIMIZE_CLASS
// clang-format off
MinimizeStyle(fire,MinFire);
// clang-format on
#else

#ifndef LMP_MIN_FIRE_H
#define LMP_MIN_FIRE_H

#include "min.h"

namespace LAMMPS_NS {

class MinFire : public Min {
 public:
  MinFire(class LAMMPS *);

  void init() override;
  void setup_style() override;
  void reset_vectors() override;
  int iterate(int) override;

 private:
  double dt, dtmax, dtmin;
  double alpha;
  bigint last_negative, ntimestep_start;
  int vdotf_negatif, flagv0;
  template <int INTEGRATOR, bool ABCFLAG> int run_iterate(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
