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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/cartesian,ComputeStressCartesian);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_CARTESIAN_H
#define LMP_COMPUTE_STRESS_CARTESIAN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressCartesian : public Compute {
 public:
  ComputeStressCartesian(class LAMMPS *, int, char **);
  ~ComputeStressCartesian() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int nbins1, nbins2, dir1, dir2, dims;
  double bin_width1, bin_width2, invV;
  bool compute_ke = true;
  bool compute_pair = true;
  bool compute_bond = true;

  // Number density, kinetic and configurational contribution to the pressure.
  double *dens, *pkxx, *pkyy, *pkzz, *pcxx, *pcyy, *pczz;
  double *tdens, *tpkxx, *tpkyy, *tpkzz, *tpcxx, *tpcyy, *tpczz;
  class NeighList *list;
  void compute_pressure(double, double, double, double, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
