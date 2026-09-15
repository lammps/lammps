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
   ILVES constraint solver (serial, exact-Newton).
   Adapted from GROMACS 2021 ILVES (LGPL-2.1), src/gromacs/mdlib/ilves.cpp and
   ilves_asym.cpp.  See ilves_graph.h for full attribution and ilves.h for the
   adaptation notes.
------------------------------------------------------------------------- */

#include "ilves.h"

#include "lammps.h"

#include <algorithm>
#include <cmath>

namespace LAMMPS_NS::ILVES {

Ilves::Ilves(LAMMPS *const _lmp, const std::vector<int> &catom1, const std::vector<int> &catom2,
             const std::vector<int> &cnode1, const std::vector<int> &cnode2,
             const std::vector<double> &cdist, const std::vector<double> &invmass) : lmp(_lmp)
{
  mol = std::make_unique<Molecule>(catom1, catom2, cnode1, cnode2, cdist, invmass);

  // reference bond vectors x_ab and predicted bond vectors xprime_ab; both are
  // needed to assemble the exact-Newton Jacobian (which uses r != s)
  for (int d = 0; d < DIM; ++d) x_ab[d].resize(mol->bonds.num);
  for (int d = 0; d < DIM; ++d) xprime_ab[d].resize(mol->bonds.num);

  // build the sparse direct solver for the constraint connectivity; it computes
  // a fill-reducing reordering, returned in solver_perm
  std::vector<int> solver_perm;
  solver = std::make_unique<SparseDirectSolver>(mol->bonds.graph, solver_perm);

  // apply the fill-reducing permutation the solver computed to the bond data
  // (the bond graph itself is no longer needed once the solver is built)
  mol->renumber_bonds(solver_perm);

  make_weights();

  // current_lagr must be as large as the rhs so the two can be swapped
  current_lagr.resize(solver->rhs.size());
}

/* ----------------------------------------------------------------------
   The fix drives each Newton iteration as prepare() then step(): prepare()
   assembles g(x) and the reference/predicted bond vectors; step() assembles and
   factors the (exact-Newton) Jacobian, solves it, and accumulates the position
   increment.  The Jacobian is re-assembled and re-factored every iteration.
------------------------------------------------------------------------- */

double Ilves::prepare(double **const x, double **const xprime)
{
  return make_rhs(x, xprime, true);
}

void Ilves::step(double **const dx)
{
  make_lhs(xprime_ab, x_ab);
  solver->LU_factor();
  solver->LU_solve();
  accumulate_increment(dx);
}

/* ----------------------------------------------------------------------
   Compute the right-hand side g(x) (the bond-length violations) into the solver
   rhs.  Returns the largest relative (squared) bond-length violation.  The bond
   vectors use raw differences of the GEOMETRY atoms (bonds.atom1/atom2): those
   are the nearest periodic image (Domain::closest_image), so the subtraction is
   already the minimum-image vector and is correct at any box size.  Graph
   connectivity uses the separate node ids (bonds.node1/node2), so a wrapped
   bond is still a single edge.
------------------------------------------------------------------------- */

double Ilves::make_rhs(double **const x, double **const xprime, const bool compute_x_ab)
{
  const int n = mol->bonds.num;
  double *const rhs = solver->rhs.data();

  double rel = 0;

  for (int k = 0; k < n; ++k) {
    const int a = mol->bonds.atom1[k];
    const int b = mol->bonds.atom2[k];

    if (compute_x_ab)
      for (int d = 0; d < DIM; ++d) x_ab[d][k] = x[b][d] - x[a][d];

    double rcd[DIM];
    for (int d = 0; d < DIM; ++d) rcd[d] = xprime[b][d] - xprime[a][d];
    for (int d = 0; d < DIM; ++d) xprime_ab[d][k] = rcd[d];

    const double scalar = rcd[XX] * rcd[XX] + rcd[YY] * rcd[YY] + rcd[ZZ] * rcd[ZZ];

    rhs[k] = 0.5 * (scalar - mol->bonds.sigma2[k]);

    rel = std::max(rel, std::abs(rhs[k]) * mol->bonds.invsigma2[k]);
  }

  return rel;
}

/* ----------------------------------------------------------------------
   Compute the left-hand side (Jacobian) from the bond vectors xab1 and xab2
   (each x_ab or xprime_ab), one stored matrix entry at a time.
------------------------------------------------------------------------- */

void Ilves::make_lhs(const BondVecs &xab1, const BondVecs &xab2)
{
  const int nentries = (int) solver->lhs.size();
  double *const lhs = solver->lhs.data();
  const int *const grows = solver->grows.data();
  const int *const gcols = solver->gcols.data();

  for (int e = 0; e < nentries; ++e) {
    const int grow = grows[e];
    const int gcol = gcols[e];

    const double scalar = xab1[XX][grow] * xab2[XX][gcol] + xab1[YY][grow] * xab2[YY][gcol] +
        xab1[ZZ][grow] * xab2[ZZ][gcol];

    lhs[e] = lhs_weights[e] * scalar;
  }
}

/* ---------------------------------------------------------------------- */

void Ilves::update_current_lagr(const bool first_time)
{
  if (first_time) {
    std::swap(current_lagr, solver->rhs);
  } else {
    const int n = mol->bonds.num;
    auto *rhs_data = solver->rhs.data();
    auto *current_lagr_data = current_lagr.data();
    for (int k = 0; k < n; ++k) current_lagr_data[k] += rhs_data[k];
  }
}

/* ----------------------------------------------------------------------
   Accumulate this Newton iteration's position increments into dx, for both
   atoms of every constraint (the atoms may be local or ghost).  dx must be
   pre-zeroed; the caller reverse-sums dx to the owning ranks and applies it.
------------------------------------------------------------------------- */

void Ilves::add_global_virial(double *const v6, const double inv_dtfsq) const
{
  // for constraint k the force on atom a is +lambda_k*inv_dtfsq*r_k (r_k = x_b -
  // x_a) and on atom b is the negative of that; the pairwise virial contribution
  // is (x_a - x_b) (x) f_a = -lambda_k*inv_dtfsq * r_k (x) r_k.
  const int n = mol->bonds.num;
  for (int k = 0; k < n; ++k) {
    const double s = -current_lagr[k] * inv_dtfsq;
    const double rx = x_ab[XX][k], ry = x_ab[YY][k], rz = x_ab[ZZ][k];
    v6[0] += s * rx * rx;
    v6[1] += s * ry * ry;
    v6[2] += s * rz * rz;
    v6[3] += s * rx * ry;
    v6[4] += s * rx * rz;
    v6[5] += s * ry * rz;
  }
}

/* ---------------------------------------------------------------------- */

double Ilves::recompute(double **const x, double **const xprime, const bool first_iter)
{
  update_current_lagr(first_iter);
  return make_rhs(x, xprime, false);
}

/* ---------------------------------------------------------------------- */

void Ilves::accumulate_increment(double **const dx) const
{
  const int n = mol->bonds.num;
  const auto *rhs = solver->rhs.data();

  for (int k = 0; k < n; ++k) {
    const int a = mol->bonds.atom1[k];
    const int b = mol->bonds.atom2[k];

    const double rhs_a = rhs[k] * mol->atoms.invmass[a];
    const double rhs_b = rhs[k] * mol->atoms.invmass[b];
    for (int d = 0; d < DIM; ++d) {
      dx[a][d] += rhs_a * x_ab[d][k];
      dx[b][d] -= rhs_b * x_ab[d][k];
    }
  }
}

/* ----------------------------------------------------------------------
   Precompute the weights of the entries of the lhs.  See the extended
   derivation in the GROMACS reference (CCKM); the weights fold in the masses
   and signs so the Jacobian A = Dg(r) inv(M) Dg(s)' can be assembled one entry
   at a time.
------------------------------------------------------------------------- */

void Ilves::make_weights()
{
  const int nentries = solver->fill_matrix.num_edges();
  lhs_weights.resize(nentries);

  for (int i = 0; i < nentries; ++i) {
    const int row = solver->grows[i];
    const int col = solver->gcols[i];

    // the shared atom of two coupled constraints is detected with the canonical
    // node ids, not the geometry indices, so an atom shared through a periodic
    // image (owner in one bond, ghost image in the other) is still recognized.
    // invmass is the same for a node and its ghost, so indexing it by node id is
    // consistent with the geometry-index use in make_rhs/accumulate_increment.
    const int arow1 = mol->bonds.node1[row];
    const int arow2 = mol->bonds.node2[row];

    if (solver->is_fillin[i]) {
      lhs_weights[i] = 0;
    } else if (row != col) {
      const int acol1 = mol->bonds.node1[col];
      const int acol2 = mol->bonds.node2[col];

      const int common = ((arow1 == acol1) || (arow1 == acol2)) ? arow1 : arow2;

      lhs_weights[i] = mol->atoms.invmass[common];

      if ((arow1 == acol2) || (arow2 == acol1)) lhs_weights[i] = -lhs_weights[i];
    } else {
      lhs_weights[i] = mol->atoms.invmass[arow1] + mol->atoms.invmass[arow2];
    }
  }
}

/* ----------------------------------------------------------------------
   estimate the solver's memory footprint in bytes
------------------------------------------------------------------------- */

double Ilves::memory_usage() const
{
  double bytes = 0.0;
  if (mol) bytes += mol->memory_usage();
  if (solver) bytes += solver->memory_usage();
  bytes += (double) (lhs_weights.size() + current_lagr.size()) * sizeof(double);
  for (int d = 0; d < DIM; ++d)
    bytes += (double) (x_ab[d].size() + xprime_ab[d].size()) * sizeof(double);
  return bytes;
}

}    // namespace LAMMPS_NS::ILVES
