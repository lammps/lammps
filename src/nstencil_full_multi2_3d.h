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

NStencilStyle(full/multi/tiered/3d,
              NStencilFullMultiTiered3d,
              NS_FULL | NS_Multi_Tiered | NS_3D |
              NS_NEWTON | NS_NEWTOFF | NS_ORTHO | NS_TRI)

#else

#ifndef LMP_NSTENCIL_FULL_MULTI_TIERED_3D_H
#define LMP_NSTENCIL_FULL_MULTI_TIERED_3D_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilFullBytype3d : public NStencil {
 public:
  NStencilFullBytype3d(class LAMMPS *);
  ~NStencilFullBytype3d();
  void create();
  void create_setup();

private:
  int ** maxstencil_type;

  int  copy_neigh_info_bytype(int);
  void create_newtoff(int, int, double);

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
