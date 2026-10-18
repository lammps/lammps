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

#include "pair_lj_cut_tip4p_long_functor_omp.h"

using namespace LAMMPS_NS;

// All behavior comes from the PairFunctorTIP4POMP<EvaluatorLJCut> threaded
// driver; this translation unit anchors the vtable and instantiates the driver
// for the lj/cut/tip4p/long/functor/omp style.

PairLJCutTIP4PLongFunctorOMP::PairLJCutTIP4PLongFunctorOMP(LAMMPS *lmp) :
    PairFunctorTIP4POMP<functor::EvaluatorLJCut>(lmp)
{
}
