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

#ifndef LMP_FIX_NH_ASPHERE_OMP_H
#define LMP_FIX_NH_ASPHERE_OMP_H

#include "fix_nh_omp.h"

namespace LAMMPS_NS {

class FixNHAsphereOMP : public FixNHOMP {
 public:
  FixNHAsphereOMP(class LAMMPS *, int, char **);
  virtual ~FixNHAsphereOMP() {}
  virtual void init();

 protected:
  double dtq;
  class AtomVecEllipsoid *avec;

  virtual void nve_v();
  virtual void nve_x();
  virtual void nh_v_temp();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Compute nvt/nph/npt asphere requires atom style ellipsoid

Self-explanatory.

E: Fix nvt/nph/npt asphere requires extended particles

The shape setting for a particle in the fix group has shape = 0.0,
which means it is a point particle.

*/
