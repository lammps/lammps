// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   Ref: Wang, Yu, Langston, Fraige, Particle shape effects in discrete
   element modelling of cohesive angular particles, Granular Matter 2011,
   13:1-12.
------------------------------------------------------------------------- */

#include "pair_body_rounded_polyhedron_omp.h"
#include "pair_body_rounded_polyhedron_omp_impl.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
template class PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedron>;
}

/* ---------------------------------------------------------------------- */

PairBodyRoundedPolyhedronOMP::PairBodyRoundedPolyhedronOMP(LAMMPS *lmp) :
    PairBodyRoundedPolyhedronOMPT<PairBodyRoundedPolyhedron>(lmp)
{
}
