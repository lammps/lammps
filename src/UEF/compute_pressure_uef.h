/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pressure/uef,ComputePressureUef);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_UEF_H
#define LMP_COMPUTE_PRESSURE_UEF_H

#include "compute_pressure.h"

namespace LAMMPS_NS {

class ComputePressureUef : public ComputePressure {
 public:
  ComputePressureUef(class LAMMPS *, int, char **);

  void init() override;
  void compute_vector() override;
  double compute_scalar() override;
  void update_rot();
  bool in_fix;    //true if this compute is used in fix/nvt/npt

 protected:
  bool ext_flags[3];    // true if used in average output pressure
  void virial_rot(double *, const double[3][3]);
  int ifix_uef;
  double rot[3][3];
};

}    // namespace LAMMPS_NS

#endif
#endif
