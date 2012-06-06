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

#ifdef PAIR_CLASS

PairStyle(hybrid/overlay,PairHybridOverlay)

#else

#ifndef LMP_PAIR_HYBRID_OVERLAY_H
#define LMP_PAIR_HYBRID_OVERLAY_H

#include "pair_hybrid.h"

namespace LAMMPS_NS {

class PairHybridOverlay : public PairHybrid {
 public:
  PairHybridOverlay(class LAMMPS *);
  ~PairHybridOverlay() {}
  void coeff(int, char **);

 private:
  void modify_requests();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair coeff for hybrid has invalid style

Style in pair coeff must have been listed in pair_style command.

*/
