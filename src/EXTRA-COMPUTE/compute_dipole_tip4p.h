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
ComputeStyle(dipole/tip4p,ComputeDipoleTIP4P);
// clang-format on
#else

#ifndef LMP_COMPUTE_DIPOLE_TIP4P_H
#define LMP_COMPUTE_DIPOLE_TIP4P_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDipoleTIP4P : public Compute {
 public:
  ComputeDipoleTIP4P(class LAMMPS *, int, char **);
  ~ComputeDipoleTIP4P() override;
  void init() override;
  void compute_vector() override;
  double compute_scalar() override;

 private:
  int usecenter;

  int typeO, typeH;
  double alpha;
  void find_M(int i, double *xM);
};

}    // namespace LAMMPS_NS

#endif
#endif
