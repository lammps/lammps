/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NSTENCIL_CLASS
// clang-format off
NStencilStyle(half/multi_old/3d,
              NStencilHalfMultiOld3d, NS_HALF | NS_MULTI_OLD | NS_3D | NS_ORTHO);
// clang-format on
#else

#ifndef LMP_NSTENCIL_HALF_MULTI_OLD_3D_H
#define LMP_NSTENCIL_HALF_MULTI_OLD_3D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfMultiOld3d : public NStencil {
 public:
  NStencilHalfMultiOld3d(class LAMMPS *);
  ~NStencilHalfMultiOld3d() {}
  void create();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
