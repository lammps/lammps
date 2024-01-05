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
typedef NStencilBinIntel<0, 0, 0> NStencilFullBin2dIntel;
NStencilStyle(full/bin/2d/intel,
              NStencilFullBin2dIntel,
              NS_FULL | NS_BIN | NS_2D | NS_ORTHO | NS_TRI | NS_INTEL);

typedef NStencilBinIntel<0, 1, 0> NStencilFullBin3dIntel;
NStencilStyle(full/bin/3d/intel,
              NStencilFullBin3dIntel,
              NS_FULL | NS_BIN | NS_3D | NS_ORTHO | NS_TRI | NS_INTEL);

typedef NStencilBinIntel<1, 0, 0> NStencilHalfBin2dIntel;
NStencilStyle(half/bin/2d/intel,
              NStencilHalfBin2dIntel,
              NS_HALF | NS_BIN | NS_2D | NS_ORTHO | NS_INTEL);

typedef NStencilBinIntel<1, 0, 1> NStencilHalfBin2dTriIntel;
NStencilStyle(half/bin/2d/tri/intel,
              NStencilHalfBin2dTriIntel,
              NS_HALF | NS_BIN | NS_2D | NS_TRI | NS_INTEL);

typedef NStencilBinIntel<1, 1, 0> NStencilHalfBin3dIntel;
NStencilStyle(half/bin/3d/intel,
              NStencilHalfBin3dIntel,
              NS_HALF | NS_BIN | NS_3D | NS_ORTHO | NS_INTEL);

typedef NStencilBinIntel<1, 1, 1> NStencilHalfBin3dTriIntel;
NStencilStyle(half/bin/3d/tri/intel,
              NStencilHalfBin3dTriIntel,
              NS_HALF | NS_BIN | NS_3D | NS_TRI | NS_INTEL);
// clang-format on
#else

#ifndef LMP_NSTENCIL_BIN_INTEL_H
#define LMP_NSTENCIL_BIN_INTEL_H

#include "nstencil.h"

namespace LAMMPS_NS {

template<int HALF, int DIM_3D, int TRI>
class NStencilBinIntel : public NStencil {
 public:
  NStencilBinIntel(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
