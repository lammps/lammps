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
typedef NStencilBin<0, 0, 0> NStencilFullBin2d;
NStencilStyle(full/bin/2d,
              NStencilFullBin2d,
              NS_FULL | NS_BIN | NS_2D | NS_ORTHO | NS_TRI);

typedef NStencilBin<0, 1, 0> NStencilFullBin3d;
NStencilStyle(full/bin/3d,
              NStencilFullBin3d,
              NS_FULL | NS_BIN | NS_3D | NS_ORTHO | NS_TRI);

typedef NStencilBin<1, 0, 0> NStencilHalfBin2d;
NStencilStyle(half/bin/2d,
              NStencilHalfBin2d,
              NS_HALF | NS_BIN | NS_2D | NS_ORTHO);

typedef NStencilBin<1, 0, 1> NStencilHalfBin2dTri;
NStencilStyle(half/bin/2d/tri,
              NStencilHalfBin2dTri,
              NS_HALF | NS_BIN | NS_2D | NS_TRI);

typedef NStencilBin<1, 1, 0> NStencilHalfBin3d;
NStencilStyle(half/bin/3d,
              NStencilHalfBin3d,
              NS_HALF | NS_BIN | NS_3D | NS_ORTHO);

typedef NStencilBin<1, 1, 1> NStencilHalfBin3dTri;
NStencilStyle(half/bin/3d/tri,
              NStencilHalfBin3dTri,
              NS_HALF | NS_BIN | NS_3D | NS_TRI);
// clang-format on
#else

#ifndef LMP_NSTENCIL_BIN_H
#define LMP_NSTENCIL_BIN_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int HALF, int DIM_3D, int TRI>
class NStencilBin : public NStencil {
 public:
  NStencilBin(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
