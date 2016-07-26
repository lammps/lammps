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

#ifdef NEIGH_STENCIL_CLASS

NeighStencilStyle(NEIGH_STENCIL_FULL_GHOST_BIN_3D,NeighStencilFullGhostBin3d)

#else

#ifndef LMP_NEIGH_STENCIL_FULL_GHOST_BIN_3D_H
#define LMP_NEIGH_STENCIL_FULL_GHOST_BIN_3D_H

#include "neigh_stencil.h"

namespace LAMMPS_NS {

class NeighStencilFullGhostBin3d : public NeighStencil {
 public:
  NeighStencilFullGhostBin3d(class LAMMPS *);
  ~NeighStencilFullGhostBin3d() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
