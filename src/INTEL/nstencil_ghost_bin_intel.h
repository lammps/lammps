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
typedef NStencilGhostBinIntel<0, 0, 0> NStencilFullGhostBin2dIntel;
NStencilStyle(full/ghost/bin/2d/intel,
              NStencilFullGhostBin2dIntel,
              NS_FULL | NS_GHOST | NS_BIN | NS_2D | NS_ORTHO | NS_TRI | NS_INTEL);

typedef NStencilGhostBinIntel<0, 1, 0> NStencilFullGhostBin3dIntel;
NStencilStyle(full/ghost/bin/3d/intel,
              NStencilFullGhostBin3dIntel,
              NS_FULL | NS_GHOST | NS_BIN | NS_3D | NS_ORTHO | NS_TRI | NS_INTEL);

typedef NStencilGhostBinIntel<1, 0, 0> NStencilHalfGhostBin2dIntel;
NStencilStyle(half/ghost/bin/2d/intel,
              NStencilHalfGhostBin2dIntel,
              NS_HALF | NS_GHOST | NS_BIN | NS_2D | NS_ORTHO | NS_INTEL);

typedef NStencilGhostBinIntel<1, 0, 1> NStencilHalfGhostBin2dTriIntel;
NStencilStyle(half/ghost/bin/2d/tri/intel,
              NStencilHalfGhostBin2dTriIntel,
              NS_HALF | NS_GHOST | NS_BIN | NS_2D | NS_TRI | NS_INTEL);

typedef NStencilGhostBinIntel<1, 1, 0> NStencilHalfGhostBin3dIntel;
NStencilStyle(half/ghost/bin/3d/intel,
              NStencilHalfGhostBin3dIntel,
              NS_HALF | NS_GHOST | NS_BIN | NS_3D | NS_ORTHO | NS_INTEL);

typedef NStencilGhostBinIntel<1, 1, 1> NStencilHalfGhostBin3dTriIntel;
NStencilStyle(half/ghost/bin/3d/tri/intel,
              NStencilHalfGhostBin3dTriIntel,
              NS_HALF | NS_GHOST | NS_BIN | NS_3D | NS_TRI | NS_INTEL);
// clang-format on
#else

#ifndef LMP_NSTENCIL_GHOST_BIN_INTEL_H
#define LMP_NSTENCIL_GHOST_BIN_INTEL_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int HALF, int DIM_3D, int TRI>
class NStencilGhostBinIntel : public NStencil {
 public:
  NStencilGhostBinIntel(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
