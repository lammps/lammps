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

NStencilStyle(full/ghost/bin/3d,
              NStencilFullGhostBin3d,
              NS_FULL | NS_GHOST | NS_BIN | NS_3D |
              NS_NEWTON | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_FULL_GHOST_BIN_3D_H
#define LMP_NSTENCIL_FULL_GHOST_BIN_3D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullGhostBin3d : public NStencil {
 public:
  NStencilFullGhostBin3d(class LAMMPS *);
  ~NStencilFullGhostBin3d() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
