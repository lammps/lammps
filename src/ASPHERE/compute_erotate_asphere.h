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
ComputeStyle(erotate/asphere,ComputeERotateAsphere);
// clang-format on
#else

#ifndef LMP_COMPUTE_EROTATE_ASPHERE_H
#define LMP_COMPUTE_EROTATE_ASPHERE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeERotateAsphere : public Compute {
 public:
  ComputeERotateAsphere(class LAMMPS *, int, char **);
  void init() override;
  double compute_scalar() override;

 private:
  double pfactor;
  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;
};

}    // namespace LAMMPS_NS

#endif
#endif
