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

NStencilStyle(full/multi2/2d,
              NStencilFullMulti22d, NS_FULL | NS_MULTI2 | NS_2D | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_FULL_MULTI2_2D_H
#define LMP_NSTENCIL_FULL_MULTI2_2D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullMulti22d : public NStencil {
 public:
  NStencilFullMulti22d(class LAMMPS *);
  ~NStencilFullMulti22d();
  void create();

 protected:
  void set_stencil_properties();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
