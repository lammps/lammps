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

#include "pair_lj_charmm_coul_long_functor_omp.h"

using namespace LAMMPS_NS;

// Behavior comes from the PairFunctorOMP<EvaluatorLJCharmm, CoulLong> threaded
// driver plus the shared PairLJCharmmCoulLongFunctorImpl 1-4 glue; this
// translation unit anchors the vtable and instantiates the threaded driver for
// the lj/charmm/coul/long/functor/omp style.

PairLJCharmmCoulLongFunctorOMP::PairLJCharmmCoulLongFunctorOMP(LAMMPS *lmp) :
    PairLJCharmmCoulLongFunctorImpl<PairFunctorOMP<functor::EvaluatorLJCharmm, functor::CoulLong>>(lmp)
{
}
