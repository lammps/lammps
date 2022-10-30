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
KSpaceStyle(msm/cg/omp,MSMCGOMP);
// clang-format on
#else

#ifndef LMP_MSM_CG_OMP_H
#define LMP_MSM_CG_OMP_H

#include "msm_omp.h"

namespace LAMMPS_NS {

class MSMCGOMP : public MSMOMP {
 public:
  MSMCGOMP(class LAMMPS *);
  ~MSMCGOMP() override;
  void settings(int, char **) override;
  void compute(int, int) override;
  double memory_usage() override;

 protected:
  int num_charged;
  int *is_charged;
  double smallq;

 protected:
  void particle_map() override;
  void make_rho() override;
  void fieldforce() override;
  void fieldforce_peratom() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
