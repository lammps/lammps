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

#include "pair_coul_long_functor.h"

using namespace LAMMPS_NS;

// All behavior is provided by the PairFunctor<EvaluatorNone, CoulLong> base
// template (a zero van der Waals term plus long-range Coulomb real space); this
// translation unit anchors the vtable and instantiates the coul/long/functor
// style.

PairCoulLongFunctor::PairCoulLongFunctor(LAMMPS *lmp) :
    PairFunctor<functor::EvaluatorNone, functor::CoulLong>(lmp)
{
}
