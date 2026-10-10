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
PairStyle(lj/charmm/coul/long/functor/omp,PairLJCharmmCoulLongFunctorOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_OMP_H
#define LMP_PAIR_LJ_CHARMM_COUL_LONG_FUNCTOR_OMP_H

#include "evaluator_lj_charmm.h"
#include "functor_coul_long.h"
#include "pair_functor_omp.h"
#include "pair_lj_charmm_coul_long_functor_impl.h"

namespace LAMMPS_NS {

// Threaded lj/charmm/coul/long/functor: the same CHARMM 1-4 glue mixin
// (PairLJCharmmCoulLongFunctorImpl) layered over the threaded driver
// PairFunctorOMP<EvaluatorLJCharmm, CoulLong> instead of the serial PairFunctor.

class PairLJCharmmCoulLongFunctorOMP :
    public PairLJCharmmCoulLongFunctorImpl<
        PairFunctorOMP<functor::EvaluatorLJCharmm, functor::CoulLong>> {
 public:
  PairLJCharmmCoulLongFunctorOMP(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
