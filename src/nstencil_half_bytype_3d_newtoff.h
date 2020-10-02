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

NStencilStyle(half/bytype/3d/newtoff,
              NStencilHalfBytype3dNewtoff,
              NS_HALF | NS_BYTYPE | NS_3D | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_HALF_BYTYPE_3D_NEWTOFF_H
#define LMP_NSTENCIL_HALF_BYTYPE_3D_NEWTOFF_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfBytype3dNewtoff : public NStencil {
 public:
  NStencilHalfBytype3dNewtoff(class LAMMPS *);
  ~NStencilHalfBytype3dNewtoff();
  void create_setup();
  void create();

 private:
  int ** maxstencil_type;

  void copy_bin_info_bytype(int);
  int  copy_neigh_info_bytype(int);
  void create_newtoff(int, int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
