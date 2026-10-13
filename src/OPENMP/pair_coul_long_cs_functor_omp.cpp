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

#include "pair_coul_long_cs_functor_omp.h"

using namespace LAMMPS_NS;

// Threaded coul/long/cs/functor; behavior comes from the
// PairFunctorOMP<EvaluatorNone, CoulLongCS> driver.  single() is disabled to
// match the serial style (the core-shell kernel does not match single).

PairCoulLongCSFunctorOMP::PairCoulLongCSFunctorOMP(LAMMPS *lmp) :
    PairFunctorOMP<functor::EvaluatorNone, functor::CoulLongCS>(lmp)
{
  single_enable = 0;
}
