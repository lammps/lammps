/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/omp,FixNVEOMP);
// clang-format on
#else

#ifndef LMP_FIX_NVE_OMP_H
#define LMP_FIX_NVE_OMP_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEOMP : public FixNVE {
 public:
  FixNVEOMP(class LAMMPS *, int, char **);

  void initial_integrate(int) override;
  void final_integrate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
