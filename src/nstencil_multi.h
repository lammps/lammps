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
typedef NStencilMulti<0, 0, 0> NStencilFullMulti2d;
NStencilStyle(full/multi/2d,
              NStencilFullMulti2d,
              NS_FULL | NS_MULTI | NS_2D | NS_ORTHO | NS_TRI);

typedef NStencilMulti<0, 1, 0> NStencilFullMulti3d;
NStencilStyle(full/multi/3d,
              NStencilFullMulti3d,
              NS_FULL | NS_MULTI | NS_3D | NS_ORTHO | NS_TRI);

typedef NStencilMulti<1, 0, 0> NStencilHalfMulti2d;
NStencilStyle(half/multi/2d,
              NStencilHalfMulti2d,
              NS_HALF | NS_MULTI | NS_2D | NS_ORTHO);

typedef NStencilMulti<1, 0, 1> NStencilHalfMulti2dTri;
NStencilStyle(half/multi/2d/tri,
              NStencilHalfMulti2dTri,
              NS_HALF | NS_MULTI | NS_2D | NS_TRI);

typedef NStencilMulti<1, 1, 0> NStencilHalfMulti3d;
NStencilStyle(half/multi/3d,
              NStencilHalfMulti3d,
              NS_HALF | NS_MULTI | NS_3D | NS_ORTHO);

typedef NStencilMulti<1, 1, 1> NStencilHalfMulti3dTri;
NStencilStyle(half/multi/3d/tri,
              NStencilHalfMulti3dTri,
              NS_HALF | NS_MULTI | NS_3D | NS_TRI);
// clang-format on
#else

#ifndef LMP_NSTENCIL_MULTI_H
#define LMP_NSTENCIL_MULTI_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int HALF, int DIM_3D, int TRI>
class NStencilMulti : public NStencil {
 public:
  NStencilMulti(class LAMMPS *);
  void create() override;

 protected:
  void set_stencil_properties() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
