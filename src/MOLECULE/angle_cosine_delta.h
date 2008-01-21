/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef ANGLE_COSINE_DELTA_H
#define ANGLE_COSINE_DELTA_H

#include "stdio.h"
#include "angle_cosine_squared.h"

namespace LAMMPS_NS {

class AngleCosineDelta : public AngleCosineSquared {
 public:
  AngleCosineDelta(class LAMMPS *);
  void compute(int, int);
  double single(int, int, int, int);
};

}

#endif
