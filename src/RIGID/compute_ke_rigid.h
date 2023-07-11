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
ComputeStyle(ke/rigid,ComputeKERigid);
// clang-format on
#else

#ifndef LMP_COMPUTE_KE_RIGID_H
#define LMP_COMPUTE_KE_RIGID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeKERigid : public Compute {
 public:
  ComputeKERigid(class LAMMPS *, int, char **);
  ~ComputeKERigid() override;
  void init() override;
  double compute_scalar() override;

 private:
  int irfix;
  char *rfix;
};

}    // namespace LAMMPS_NS

#endif
#endif
