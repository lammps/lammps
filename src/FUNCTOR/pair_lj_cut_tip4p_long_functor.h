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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/tip4p/long/functor,PairLJCutTIP4PLongFunctor);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_FUNCTOR_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_FUNCTOR_H

#include "evaluator_lj_cut.h"
#include "pair_functor_tip4p.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongFunctor : public PairFunctorTIP4P<functor::EvaluatorLJCut> {
 public:
  PairLJCutTIP4PLongFunctor(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
