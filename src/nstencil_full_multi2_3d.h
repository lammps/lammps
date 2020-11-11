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

NStencilStyle(full/multi2/3d,
              NStencilFullMulti23d, NS_FULL | NS_Multi2 | NS_3D | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_FULL_MULTI2_3D_H
#define LMP_NSTENCIL_FULL_MULTI2_3D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullMulti23d : public NStencil {
 public:
  NStencilFullMulti23d(class LAMMPS *);
  ~NStencilFullMulti23d();
  void create();

 protected:
  void set_stencil_properties();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
