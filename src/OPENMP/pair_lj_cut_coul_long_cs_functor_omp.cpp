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

#include "pair_lj_cut_coul_long_cs_functor_omp.h"

using namespace LAMMPS_NS;

// Threaded lj/cut/coul/long/cs/functor; all behavior comes from the
// PairFunctorOMP<EvaluatorLJCut, CoulLongCS> driver.

PairLJCutCoulLongCSFunctorOMP::PairLJCutCoulLongCSFunctorOMP(LAMMPS *lmp) :
    PairFunctorOMP<functor::EvaluatorLJCut, functor::CoulLongCS>(lmp)
{
  // single() uses the inherited A-polynomial Coulomb and does not match the
  // B-polynomial compute kernel (as in the original lj/cut/coul/long/cs)
  single_enable = 0;
}
