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

NStencilStyle(cac,
              NStencilCAC,
              NS_CAC | NS_BIN | NS_3D |
              NS_NEWTON | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_CAC_H
#define LMP_NSTENCIL_CAC_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilCAC : public NStencil {
 public:
  NStencilCAC(class LAMMPS *);
  ~NStencilCAC() {}
  void create();
  void create_setup();
  void post_create();
  void post_create_setup();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
