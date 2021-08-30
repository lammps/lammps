/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NSTENCIL_SSA_H
#define LMP_NSTENCIL_SSA_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilSSA : public NStencil {
 public:
  NStencilSSA(class LAMMPS *lmp) : NStencil(lmp) { xyzflag = 1; }
  ~NStencilSSA() {}
  virtual void create() = 0;

  // first stencil index for each subphase, with last index at end
  int nstencil_ssa[5];
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

*/
