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

#ifndef LMP_FIX_NH_SPHERE_OMP_H
#define LMP_FIX_NH_SPHERE_OMP_H

#include "fix_nh_omp.h"

namespace LAMMPS_NS {

class FixNHSphereOMP : public FixNHOMP {
 public:
  FixNHSphereOMP(class LAMMPS *, int, char **);
  virtual ~FixNHSphereOMP() {}
  virtual void init();

 protected:
  virtual void nve_v();
  virtual void nh_v_temp();
};

}

#endif

/* ERROR/WARNING messages:

E: Fix nvt/nph/npt sphere requires atom style sphere

Self-explanatory.

E: Fix nvt/sphere requires extended particles

This fix can only be used for particles of a finite size.

*/
