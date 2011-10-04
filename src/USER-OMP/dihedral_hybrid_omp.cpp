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

#include "dihedral_hybrid_omp.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DihedralHybridOMP::DihedralHybridOMP(LAMMPS *lmp) : DihedralHybrid(lmp)
{
  suffix = "omp";
}

