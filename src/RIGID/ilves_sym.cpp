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
   ILVES-FAST (symmetric variant).  Adapted from GROMACS 2021 ILVES (LGPL-2.1),
   src/gromacs/mdlib/ilves_sym.cpp.  The serial single-partition solve drops the
   OpenMP/MPI/SIMD machinery of the reference; see ilves_graph.h for attribution.
------------------------------------------------------------------------- */

#include "ilves_sym.h"

#include <cmath>
#include <utility>

namespace LAMMPS_NS {
namespace ILVES {

IlvesSym::IlvesSym(LAMMPS *const lmp, const int nbonds, const int *const catom1,
                   const int *const catom2, const real *const cdist, const real *const invmass,
                   const int nthreads) :
    Ilves(lmp, nbonds, catom1, catom2, cdist, invmass, nthreads, true)
{
}

std::pair<bool, int> IlvesSym::solve(double **const x, double **const xprime, const real tol,
                                     const int maxiter)
{
  int numit = 0;

  // Compute g(x); also fills the reference bond vectors x_ab.
  real ptau = make_rhs(0, x, xprime, true);

  // Assemble the (symmetric, step-constant) Jacobian and factor it once.
  make_lhs(0, x_ab, x_ab);
  schur_solver->LDLT_factor();

  // Newton iteration: reuse the cached factor across iterations.
  for (int i = 0; (i < maxiter) && std::isfinite(ptau) && (tol < ptau); ++i) {
    ++numit;

    schur_solver->LDLT_solve();
    update_positions(0, xprime);
    update_current_lagr(0, i == 0);

    ptau = make_rhs(0, x, xprime, false);
  }

  return std::make_pair(std::isfinite(ptau) && (ptau < tol), numit);
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
