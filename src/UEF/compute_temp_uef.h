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
ComputeStyle(temp/uef,ComputeTempUef);
// clang-format on
#else

#ifndef LMP_COMPUTE_TEMP_UEF_H
#define LMP_COMPUTE_TEMP_UEF_H

#include "compute_temp.h"

namespace LAMMPS_NS {

class ComputeTempUef : public ComputeTemp {
 public:
  ComputeTempUef(class LAMMPS *, int, char **);

  void init() override;
  void compute_vector() override;
  void yes_rot();
  void no_rot();

 protected:
  bool rot_flag;
  void virial_rot(double *, const double[3][3]);
  int ifix_uef;
};

}    // namespace LAMMPS_NS

#endif
#endif
