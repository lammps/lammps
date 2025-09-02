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

#ifdef NSTENCIL_CLASS
// clang-format off
using NStencilFullMultiOld2d = NStencilMultiOld<0, 0, 0>;
NStencilStyle(full/multi/old/2d,
              NStencilFullMultiOld2d,
              NS_FULL | NS_MULTI_OLD | NS_2D | NS_ORTHO | NS_TRI);

using NStencilFullMultiOld3d = NStencilMultiOld<0, 1, 0>;
NStencilStyle(full/multi/old/3d,
              NStencilFullMultiOld3d,
              NS_FULL | NS_MULTI_OLD | NS_3D | NS_ORTHO | NS_TRI);

using NStencilHalfMultiOld2d = NStencilMultiOld<1, 0, 0>;
NStencilStyle(half/multi/old/2d,
              NStencilHalfMultiOld2d,
              NS_HALF | NS_MULTI_OLD | NS_2D | NS_ORTHO);

using NStencilHalfMultiOld2dTri = NStencilMultiOld<1, 0, 1>;
NStencilStyle(half/multi/old/2d/tri,
              NStencilHalfMultiOld2dTri,
              NS_HALF | NS_MULTI_OLD | NS_2D | NS_TRI);

using NStencilHalfMultiOld3d = NStencilMultiOld<1, 1, 0>;
NStencilStyle(half/multi/old/3d,
              NStencilHalfMultiOld3d,
              NS_HALF | NS_MULTI_OLD | NS_3D | NS_ORTHO);

using NStencilHalfMultiOld3dTri = NStencilMultiOld<1, 1, 1>;
NStencilStyle(half/multi/old/3d/tri,
              NStencilHalfMultiOld3dTri,
              NS_HALF | NS_MULTI_OLD | NS_3D | NS_TRI);
// clang-format on
#else

#ifndef LMP_NSTENCIL_MULTI_OLD_H
#define LMP_NSTENCIL_MULTI_OLD_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int HALF, int DIM_3D, int TRI>
class NStencilMultiOld : public NStencil {
 public:
  NStencilMultiOld(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
