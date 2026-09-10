/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(orient/eco/omp,FixOrientECOOMP);
// clang-format on
#else

#ifndef LMP_FIX_ORIENT_ECO_OMP_H
#define LMP_FIX_ORIENT_ECO_OMP_H

#include "fix_orient_eco.h"

namespace LAMMPS_NS {

class FixOrientECOOMP : public FixOrientECO {
 public:
  FixOrientECOOMP(class LAMMPS *, int, char **);
  void post_force(int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
