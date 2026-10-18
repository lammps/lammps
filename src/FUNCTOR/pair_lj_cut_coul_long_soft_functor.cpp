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

#include "pair_lj_cut_coul_long_soft_functor.h"

using namespace LAMMPS_NS;

// All behavior is provided by the PairFunctor<EvaluatorLJCutSoft, CoulLongSoft>
// base template (soft-core LJ plus soft-core long-range Coulomb real space);
// this translation unit anchors the vtable and instantiates the
// lj/cut/coul/long/soft/functor style.

PairLJCutCoulLongSoftFunctor::PairLJCutCoulLongSoftFunctor(LAMMPS *lmp) :
    PairFunctor<functor::EvaluatorLJCutSoft, functor::CoulLongSoft>(lmp)
{
}
