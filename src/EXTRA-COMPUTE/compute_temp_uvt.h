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
ComputeStyle(temp/uvt,ComputeTempUVT);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_UVT_H
#define LMP_COMPUTE_TEMP_UVT_H

#include "compute_temp.h"

namespace LAMMPS_NS {

class ComputeTempUVT : public ComputeTemp {
 public:
  ComputeTempUVT(class LAMMPS *, int, char **);
  ~ComputeTempUVT() override;
  void init() override;
  double compute_scalar() override;

 protected:
  void dof_compute() override;

  char *id_fix;
  double *ne_dot, *ne_mass;
};

}    // namespace LAMMPS_NS

#endif
#endif
