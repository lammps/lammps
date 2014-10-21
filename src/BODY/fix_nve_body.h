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

#ifdef FIX_CLASS

FixStyle(nve/body,FixNVEBody)

#else

#ifndef LMP_FIX_NVE_BODY_H
#define LMP_FIX_NVE_BODY_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEBody : public FixNVE {
 public:
  FixNVEBody(class LAMMPS *, int, char **);
  void init();
  void initial_integrate(int);
  void final_integrate();

 private:
  double dtq;
  class AtomVecBody *avec;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Fix nve/body requires atom style body

Self-explanatory.

E: Fix nve/body requires bodies

This fix can only be used for particles that are bodies.

*/
