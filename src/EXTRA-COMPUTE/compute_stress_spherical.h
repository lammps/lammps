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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/spherical,ComputeStressSpherical);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_SPHERICAL_H
#define LMP_COMPUTE_STRESS_SPHERICAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressSpherical : public Compute {
 public:
  ComputeStressSpherical(class LAMMPS *, int, char **);
  ~ComputeStressSpherical() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int nbins;
  double bin_width, x0, y0, z0, Rmax;

  // Number density, kinetic and configurational contribution to the pressure.
  double *invV, *dens, *pkrr, *pktt, *pkpp, *pcrr, *pctt, *pcpp;
  double *tdens, *tpkrr, *tpktt, *tpkpp, *tpcrr, *tpctt, *tpcpp;
  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
