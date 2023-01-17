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
ComputeStyle(PRESSURE/GREM,ComputePressureGrem);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_GREM_H
#define LMP_COMPUTE_PRESSURE_GREM_H

#include "compute_pressure.h"

namespace LAMMPS_NS {

class ComputePressureGrem : public ComputePressure {
 public:
  ComputePressureGrem(class LAMMPS *, int, char **);
  ~ComputePressureGrem() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;

 protected:
  // Access to gREM fix scale factor
  char *fix_grem;
  double *scale_grem;
};

}    // namespace LAMMPS_NS

#endif
#endif
