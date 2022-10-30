/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(brownian,FixBrownian);
// clang-format on
#else

#ifndef LMP_FIX_BROWNIAN_H
#define LMP_FIX_BROWNIAN_H

#include "fix_brownian_base.h"

namespace LAMMPS_NS {

class FixBrownian : public FixBrownianBase {
 public:
  FixBrownian(class LAMMPS *, int, char **);

  void init() override;
  void initial_integrate(int) override;

 private:
  template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D> void initial_integrate_templated();
};

}    // namespace LAMMPS_NS
#endif
#endif
