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

#ifndef LMP_FIX_NH_BODY_H
#define LMP_FIX_NH_BODY_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNHBody : public FixNH {
 public:
  FixNHBody(class LAMMPS *, int, char **);
  virtual ~FixNHBody() {}
  void init();

 protected:
  double dtq;
  class AtomVecBody *avec;

  void nve_v();
  void nve_x();
  void nh_v_temp();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Compute nvt/nph/npt body requires atom style body

Self-explanatory.

E: Fix nvt/nph/npt body requires bodies

Self-explanatory.

*/
