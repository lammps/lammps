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

/* ----------------------------------------------------------------------
   ILVES (full / exact-Newton variant).  Adapted from GROMACS 2021 ILVES
   (LGPL-2.1), src/gromacs/mdlib/ilves_asym.cpp.  The serial single-partition
   solve drops the OpenMP/MPI/SIMD machinery; see ilves_graph.h for attribution.
------------------------------------------------------------------------- */

#include "ilves_asym.h"

namespace LAMMPS_NS {
namespace ILVES {

IlvesAsym::IlvesAsym(LAMMPS *const lmp, const int nbonds, const int *const catom1,
                     const int *const catom2, const double *const cdist, const double *const invmass,
                     const int nthreads) :
    Ilves(lmp, nbonds, catom1, catom2, cdist, invmass, nthreads, false)
{
  // The asymmetric variant also needs the predicted bond vectors xprime_ab.
  for (int d = 0; d < DIM; ++d) xprime_ab[d].resize(mol->bonds.num);
}

double IlvesAsym::prepare(double **const x, double **const xprime)
{
  // Compute g(x); fills the reference bond vectors x_ab and (since xprime_ab is
  // allocated for this variant) the predicted bond vectors xprime_ab.  The
  // exact-Newton Jacobian is reassembled and refactored each step (in step()).
  return make_rhs(0, x, xprime, true);
}

void IlvesAsym::step(double **const dx)
{
  make_lhs(0, xprime_ab, x_ab);
  schur_solver->LU_factor();
  schur_solver->LU_solve();
  accumulate_increment(0, dx);
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
