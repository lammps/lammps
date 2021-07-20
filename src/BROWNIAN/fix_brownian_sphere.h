/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(brownian/sphere,FixBrownianSphere);
// clang-format on
#else

#ifndef LMP_FIX_BROWNIAN_SPHERE_H
#define LMP_FIX_BROWNIAN_SPHERE_H

#include "fix_brownian_base.h"

namespace LAMMPS_NS {

class FixBrownianSphere : public FixBrownianBase {
 public:
  FixBrownianSphere(class LAMMPS *, int, char **);
  virtual ~FixBrownianSphere(){};
  void init();
  void initial_integrate(int);

 private:
  template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D> void initial_integrate_templated();
  double g3, g4;
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix brownian/sphere command.

Wrong number/type of input arguments.

E: Compute brownian/sphere requires atom style sphere

Self-explanatory.

E: Compute brownian/sphere requires atom attribute mu

Self-explanatory.

E: Fix brownian/sphere translational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix brownian/sphere rotational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix brownian/sphere translational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix brownian/sphere rotational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix brownian/sphere seed must be > 0.

*/
