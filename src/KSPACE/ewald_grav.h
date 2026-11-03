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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(ewald/grav,EwaldGrav);
// clang-format on
#else

#ifndef LMP_EWALD_GRAV_H
#define LMP_EWALD_GRAV_H

#include "ewald.h"

namespace LAMMPS_NS {

class EwaldGrav : public Ewald {
 public:
  EwaldGrav(class LAMMPS *);
  ~EwaldGrav() override;
  void init() override;
  void compute(int, int) override;
  void qsum_qsq(int warning_flag = 1) override;
  void two_charge() override;
  double memory_usage() override;

 protected:
  double gconst;     // gravitational constant G, taken from the pair style
  double *mass_src;  // per-atom masses of owned atoms
  int nmax_src;

  double *source_array() override { return mass_src; }
  void update_masses();
};

}    // namespace LAMMPS_NS

#endif
#endif
