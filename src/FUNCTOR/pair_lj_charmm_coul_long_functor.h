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
PairStyle(lj/charmm/coul/long/functor,PairLJCharmmCoulLongFunctor);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_H
#define LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_H

#include "evaluator_lj_charmm.h"
#include "functor_coul_long.h"
#include "pair_functor.h"
#include "pair_lj_charmm_coul_long_functor_impl.h"

namespace LAMMPS_NS {

// The pairwise kernel comes from PairFunctor<EvaluatorLJCharmm, CoulLong>; the
// CHARMM 1-4 (dihedral) glue is supplied by the shared
// PairLJCharmmCoulLongFunctorImpl mixin (also reused by the /omp variant).

class PairLJCharmmCoulLongFunctor :
    public PairLJCharmmCoulLongFunctorImpl<
        PairFunctor<functor::EvaluatorLJCharmm, functor::CoulLong>> {
 public:
  PairLJCharmmCoulLongFunctor(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
