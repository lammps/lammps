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
MinimizeStyle(spin,MinSpin);
// clang-format on
#else

#ifndef LMP_MIN_SPIN_H
#define LMP_MIN_SPIN_H

#include "min.h"

namespace LAMMPS_NS {

class MinSpin : public Min {
 public:
  MinSpin(class LAMMPS *);

  void init() override;
  void setup_style() override;
  int modify_param(int, char **) override;
  void reset_vectors() override;
  int iterate(int) override;
  double evaluate_dt();
  void advance_spins(double);

 private:
  // global and spin timesteps

  double dt;
  double dts;

  double alpha_damp;         // damping for spin minimization
  double discrete_factor;    // factor for spin timestep evaluation

  double *spvec;    // variables for atomic dof, as 1d vector
  double *fmvec;    // variables for atomic dof, as 1d vector

  bigint last_negative;
};

}    // namespace LAMMPS_NS

#endif
#endif
