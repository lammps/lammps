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

namespace LAMMPS_NS {
namespace ILVES {

IlvesSym::IlvesSym(LAMMPS *const lmp, const int nbonds, const int *const catom1,
                   const int *const catom2, const real *const cdist, const real *const invmass,
                   const int nthreads) :
    Ilves(lmp, nbonds, catom1, catom2, cdist, invmass, nthreads, true)
{
}

real IlvesSym::prepare(double **const x, double **const xprime)
{
  // Compute g(x) (also fills the reference bond vectors x_ab), then assemble the
  // symmetric, step-constant Jacobian and factor it once (reused every step).
  const real ptau = make_rhs(0, x, xprime, true);
  make_lhs(0, x_ab, x_ab);
  schur_solver->LDLT_factor();
  return ptau;
}

void IlvesSym::step(double **const dx)
{
  schur_solver->LDLT_solve();
  accumulate_increment(0, dx);
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
