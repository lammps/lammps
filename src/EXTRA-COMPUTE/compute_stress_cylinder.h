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
ComputeStyle(stress/cylinder,ComputeStressCylinder);
ComputeStyle(pressure/cylinder,ComputeStressCylinder);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_CYLINDER_H
#define LMP_COMPUTE_STRESS_CYLINDER_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressCylinder : public Compute {
 public:
  ComputeStressCylinder(class LAMMPS *, int, char **);
  ~ComputeStressCylinder() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int kinetic_flag;
  int nbins, nphi, nzbins;
  double *Pvr_temp, *Pvr_all, *Pvz_temp, *Pvz_all, *Pvphi_temp, *Pvphi_all;
  double *Pkr_temp, *Pkr_all, *Pkz_temp, *Pkz_all, *Pkphi_temp, *Pkphi_all;
  double *R, *Rinv, *R2, *PrAinv, *PzAinv, PphiAinv;
  double Rmax, bin_width, nktv2p;
  double *R2kin, *density_temp, *invVbin, *density_all;
  double *tangent, *ephi_x, *ephi_y;
  double *binz;

  double zlo, zhi;

  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
