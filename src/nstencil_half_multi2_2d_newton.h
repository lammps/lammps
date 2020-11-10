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

NStencilStyle(half/bytype/2d/newton,
              NStencilHalfBytype2dNewton,
              NS_HALF | NS_BYTYPE | NS_2D | NS_NEWTON | NS_ORTHO)

#else

#ifndef LMP_NSTENCIL_HALF_BYTYPE_2D_NEWTON_H
#define LMP_NSTENCIL_HALF_BYTYPE_2D_NEWTON_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilHalfBytype2dNewton : public NStencil {
 public:
  NStencilHalfBytype2dNewton(class LAMMPS *);
  ~NStencilHalfBytype2dNewton();
  void create_setup();
  void create();

 private:
  int ** maxstencil_type;

  void copy_bin_info_bytype(int);
  int  copy_neigh_info_bytype(int);
  void create_newton(int, int, double);
  void create_newtoff(int, int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
