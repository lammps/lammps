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
FixStyle(langevin/spin,FixLangevinSpin);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_SPIN_H
#define LMP_FIX_LANGEVIN_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinSpin : public Fix {
 public:
  int tdamp_flag, temp_flag;    // damping and temperature flags

  FixLangevinSpin(class LAMMPS *, int, char **);
  ~FixLangevinSpin() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void add_tdamping(double *, double *);    // add transverse damping
  void add_temperature(double[3]);
  void compute_single_langevin(int, double *, double *);

 protected:
  double alpha_t;       // transverse mag. damping
  double dts;           // magnetic timestep
  double temp;          // spin bath temperature
  double D, sigma;      // bath intensity var.
  double gil_factor;    // gilbert's prefactor

  int nlevels_respa;
  class RanMars *random;
  int seed;
};

}    // namespace LAMMPS_NS

#endif
#endif
