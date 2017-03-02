/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NSTENCIL_CLASS

NStencilStyle(half/bin/3d/newton/ssa,
              NStencilHalfBin3dNewtonSSA,
              NS_HALF | NS_BIN | NS_3D | NS_NEWTON | NS_SSA | NS_ORTHO | NS_GHOST)

#else

#ifndef LMP_NSTENCIL_HALF_BIN_3D_NEWTON_SSA_H
#define LMP_NSTENCIL_HALF_BIN_3D_NEWTON_SSA_H

#include "nstencil_ssa.h"

namespace LAMMPS_NS {

class NStencilHalfBin3dNewtonSSA : public NStencilSSA {
 public:
  NStencilHalfBin3dNewtonSSA(class LAMMPS *);
  ~NStencilHalfBin3dNewtonSSA() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
