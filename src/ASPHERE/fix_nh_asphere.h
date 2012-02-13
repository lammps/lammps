/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_ASPHERE_H
#define LMP_FIX_NH_ASPHERE_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHAsphere : public FixNH {
 public:
  FixNHAsphere(class LAMMPS *, int, char **);
  virtual ~FixNHAsphere() {}
  void init();

 protected:
  double dtq;
  class AtomVecEllipsoid *avec;

  void nve_v();
  void nve_x();
  void nh_v_temp();
};

}

#endif

/* ERROR/WARNING messages:

E: Compute nvt/nph/npt asphere requires atom style ellipsoid

Self-explanatory.

E: Fix nvt/nph/npt asphere requires extended particles

The shape setting for a particle in the fix group has shape = 0.0,
which means it is a point particle.

*/
