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
PairStyle(lj/cut/tip4p/long/functor/omp,PairLJCutTIP4PLongFunctorOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_FUNCTOR_OMP_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_FUNCTOR_OMP_H

#include "evaluator_lj_cut.h"
#include "pair_functor_tip4p_omp.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongFunctorOMP : public PairFunctorTIP4POMP<functor::EvaluatorLJCut> {
 public:
  PairLJCutTIP4PLongFunctorOMP(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
