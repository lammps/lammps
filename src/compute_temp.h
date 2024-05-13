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
ComputeStyle(temp,ComputeTemp);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_H
#define LMP_COMPUTE_TEMP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTemp : public Compute {
 public:
  ComputeTemp(class LAMMPS *, int, char **);
  ~ComputeTemp() override;
  void init() override {}
  void setup() override;
  double compute_scalar() override;
  void compute_vector() override;

 protected:
  double tfactor;

  virtual void dof_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif
