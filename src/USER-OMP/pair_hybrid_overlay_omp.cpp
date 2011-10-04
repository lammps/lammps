/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_hybrid_overlay_omp.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridOverlayOMP::PairHybridOverlayOMP(LAMMPS *lmp) : PairHybridOverlay(lmp)
{
  suffix = "omp";
}

