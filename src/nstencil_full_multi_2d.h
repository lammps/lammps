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

NStencilStyle(full/multi/2d,
              NStencilFullMulti2d,
              NS_FULL | NS_MULTI | NS_2D |
              NS_NEWTON | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_FULL_MULTI_2D_H
#define LMP_NSTENCIL_FULL_MULTI_2D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullMulti2d : public NStencil {
 public:
  NStencilFullMulti2d(class LAMMPS *);
  ~NStencilFullMulti2d() {}
  void create();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
