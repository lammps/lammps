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

#ifdef FIX_CLASS

FixStyle(nph/asphere,FixNPHAsphere)

#else

#ifndef LMP_FIX_NPH_ASPHERE_H
#define LMP_FIX_NPH_ASPHERE_H

#include "fix_nh_asphere.h"

namespace LAMMPS_NS {

class FixNPHAsphere : public FixNHAsphere {
 public:
  FixNPHAsphere(class LAMMPS *, int, char **);
  ~FixNPHAsphere() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control can not be used with fix nph/asphere

Self-explanatory.

E: Pressure control must be used with fix nph/asphere

Self-explanatory.

*/
