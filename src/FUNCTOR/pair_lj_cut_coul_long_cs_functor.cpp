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

#include "pair_lj_cut_coul_long_cs_functor.h"

using namespace LAMMPS_NS;

// All behavior is provided by the PairFunctor<EvaluatorLJCut, CoulLongCS> base
// template (the Lennard-Jones evaluator plus the core-shell long-range Coulomb
// kernel); this translation unit anchors the vtable and instantiates the
// lj/cut/coul/long/cs/functor style.

PairLJCutCoulLongCSFunctor::PairLJCutCoulLongCSFunctor(LAMMPS *lmp) :
    PairFunctor<functor::EvaluatorLJCut, functor::CoulLongCS>(lmp)
{
  // as in the original lj/cut/coul/long/cs: single() uses the inherited
  // A-polynomial Coulomb and so does not match the B-polynomial compute kernel
  single_enable = 0;
}
