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

NStencilStyle(half/multi2/2d,
              NStencilHalfMulti22d, NS_HALF | NS_MULTI2 | NS_2D | NS_ORTHO)

#else

#ifndef LMP_NSTENCIL_HALF_MULTI2_2D_H
#define LMP_NSTENCIL_HALF_MULTI2_2D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfMulti22d : public NStencil {
 public:
  NStencilHalfMulti22d(class LAMMPS *);
  ~NStencilHalfMulti22d();
  void create();

 protected:
  void set_stencil_properties();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
