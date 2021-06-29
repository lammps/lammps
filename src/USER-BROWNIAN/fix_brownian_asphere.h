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
FixStyle(brownian/asphere,FixBrownianAsphere);
// clang-format on
#else

#ifndef LMP_FIX_BROWNIAN_ASPHERE_H
#define LMP_FIX_BROWNIAN_ASPHERE_H

#include "fix_brownian_base.h"

namespace LAMMPS_NS {

class FixBrownianAsphere : public FixBrownianBase {
 public:
  FixBrownianAsphere(class LAMMPS *, int, char **);
  virtual ~FixBrownianAsphere(){};
  void initial_integrate(int);
  void init();

 protected:
  class AtomVecEllipsoid *avec;

 private:
  template <int Tp_UNIFORM, int Tp_GAUSS, int Tp_DIPOLE, int Tp_2D>
  void initial_integrate_templated();
};
}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix brownian/asphere command.

Wrong number/type of input arguments.

E: Compute brownian/asphere requires atom style sphere

Self-explanatory.

E: Compute brownian/asphere requires atom style ellipsoid

Self-explanatory.

E: Compute brownian/asphere dipole requires atom attribute mu

Self-explanatory.

E: Fix brownian/asphere translational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix brownian/asphere rotational viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix brownian/asphere translational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix brownian/asphere rotational diffusion coefficient must be > 0.

Self-explanatory.

E: Fix brownian/asphere seed must be > 0.

*/
