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
typedef NStencilGhostBin<0> NStencilFullGhostBin2d;
NStencilStyle(full/ghost/bin/2d,
              NStencilFullGhostBin2d,
              NS_FULL | NS_GHOST | NS_BIN | NS_2D | NS_ORTHO | NS_TRI);

typedef NStencilGhostBin<1> NStencilFullGhostBin3d;
NStencilStyle(full/ghost/bin/3d,
              NStencilFullGhostBin3d,
              NS_FULL | NS_GHOST | NS_BIN | NS_3D | NS_ORTHO | NS_TRI);
// clang-format on
#else

#ifndef LMP_NSTENCIL_GHOST_BIN_H
#define LMP_NSTENCIL_GHOST_BIN_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int DIM_3D>
class NStencilGhostBin : public NStencil {
 public:
  NStencilGhostBin(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
