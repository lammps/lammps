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

NStencilStyle(half/bin/2d/newton/tri,
              NStencilHalfBin2dNewtonTri,
              NS_HALF | NS_BIN | NS_2D | NS_NEWTON | NS_TRI)

#else

#ifndef LMP_NSTENCIL_HALF_BIN_2D_NEWTON_TRI_H
#define LMP_NSTENCIL_HALF_BIN_2D_NEWTON_TRI_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfBin2dNewtonTri : public NStencil {
 public:
  NStencilHalfBin2dNewtonTri(class LAMMPS *);
  ~NStencilHalfBin2dNewtonTri() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
