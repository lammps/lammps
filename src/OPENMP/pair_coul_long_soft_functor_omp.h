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
PairStyle(coul/long/soft/functor/omp,PairCoulLongSoftFunctorOMP);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_LONG_SOFT_FUNCTOR_OMP_H
#define LMP_PAIR_COUL_LONG_SOFT_FUNCTOR_OMP_H

#include "evaluator_none_soft.h"
#include "functor_coul_soft.h"
#include "pair_functor_omp.h"

namespace LAMMPS_NS {

class PairCoulLongSoftFunctorOMP :
    public PairFunctorOMP<functor::EvaluatorNoneSoft, functor::CoulLongSoft> {
 public:
  PairCoulLongSoftFunctorOMP(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
