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
NStencilStyle(half/bin/2d/ssa,
              NStencilHalfBin2dSSA,
              NS_HALF | NS_BIN | NS_2D | NS_SSA | NS_ORTHO | NS_GHOST);
// clang-format on
#else

#ifndef LMP_NSTENCIL_HALF_BIN_2D_SSA_H
#define LMP_NSTENCIL_HALF_BIN_2D_SSA_H

#include "nstencil_ssa.h"

namespace LAMMPS_NS {

class NStencilHalfBin2dSSA : public NStencilSSA {
 public:
  NStencilHalfBin2dSSA(class LAMMPS *);
  void create() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
