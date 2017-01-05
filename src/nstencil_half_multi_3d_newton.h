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

NStencilStyle(half/multi/3d/newton,
              NStencilHalfMulti3dNewton,
              NS_HALF | NS_MULTI | NS_3D | NS_NEWTON | NS_ORTHO)

#else

#ifndef LMP_NSTENCIL_HALF_MULTI_3D_NEWTON_H
#define LMP_NSTENCIL_HALF_MULTI_3D_NEWTON_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfMulti3dNewton : public NStencil {
 public:
  NStencilHalfMulti3dNewton(class LAMMPS *);
  ~NStencilHalfMulti3dNewton() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
