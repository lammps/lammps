// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_charmm_coul_long_functor.h"

using namespace LAMMPS_NS;

// All behavior is provided by the PairFunctor<EvaluatorLJCharmm, CoulLong> base
// template; this translation unit anchors the vtable and instantiates the
// templated driver for the lj/charmm/coul/long/functor style.

PairLJCharmmCoulLongFunctor::PairLJCharmmCoulLongFunctor(LAMMPS *lmp) :
    PairFunctor<functor::EvaluatorLJCharmm, functor::CoulLong>(lmp)
{
  implicit = 0;
  lj14_1 = lj14_2 = lj14_3 = lj14_4 = nullptr;
  mix_flag = ARITHMETIC;    // CHARMM uses arithmetic (Lorentz-Berthelot) mixing by default
}
