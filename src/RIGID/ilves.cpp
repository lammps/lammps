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
   ILVES constraint solver: abstract base class.
   Adapted from GROMACS 2021 ILVES (LGPL-2.1), src/gromacs/mdlib/ilves.cpp.
   See ilves_graph.h for full attribution and ilves.h for the adaptation notes.
------------------------------------------------------------------------- */

#include "ilves.h"

#include "lammps.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace LAMMPS_NS {
namespace ILVES {

Ilves::Ilves(LAMMPS *const _lmp, const int nbonds, const int *const catom1, const int *const catom2,
             const real *const cdist, const real *const invmass, const int threads,
             const bool upper_tri) :
    lmp(_lmp), nthreads(threads)
{
  mol = std::unique_ptr<Molecule>(new Molecule(nbonds, catom1, catom2, cdist, invmass));

  for (int d = 0; d < DIM; ++d) x_ab[d].resize(mol->bonds.num);

  const int nparts = nthreads + 1;

  std::vector<int> partitioning(nparts + 1);
  std::vector<int> bonds_perm_kpart;

  // Partitioning of the bonds between threads.
  if (nthreads > 1) {
    std::vector<int> bonds_part_id;
    std::vector<int> nbonds_per_part(nparts, 0);

    const int block_max_size = std::min(10, mol->bonds.num / nthreads);
    const bool is_disjoint_mol = disjoint_mol(block_max_size);

    if (is_disjoint_mol) {
      bonds_part_id = mol->bonds.graph.kway_partition_disjoint(nthreads);
      for (int bond = 0; bond < mol->bonds.num; ++bond) ++nbonds_per_part[bonds_part_id[bond]];
    } else {
      const auto atoms_part_id = mol->atoms.graph.kway_partition(nthreads);
      bonds_part_id.resize(mol->bonds.num);
      for (int bond = 0; bond < mol->bonds.graph.num_nodes(); ++bond) {
        const int la = mol->bonds.latom1[bond];
        const int lb = mol->bonds.latom2[bond];
        const int part = (atoms_part_id[la] != atoms_part_id[lb]) ? nthreads : atoms_part_id[la];
        bonds_part_id[bond] = part;
        ++nbonds_per_part[part];
      }
    }

    partitioning[0] = 0;
    for (int part = 0; part < nparts; ++part)
      partitioning[part + 1] = partitioning[part] + nbonds_per_part[part];

    std::vector<int> bonds_part_idx(nparts, 0);
    for (int p = 1; p < nparts; ++p)
      bonds_part_idx[p] = bonds_part_idx[p - 1] + nbonds_per_part[p - 1];

    bonds_perm_kpart.resize(mol->bonds.num);
    std::vector<int> bonds_iperm_kpart(mol->bonds.num);
    for (int bond = 0; bond < mol->bonds.num; ++bond) {
      const int part = bonds_part_id[bond];
      bonds_iperm_kpart[bond] = bonds_part_idx[part];
      ++bonds_part_idx[part];
    }
    for (int bond = 0; bond < mol->bonds.num; ++bond)
      bonds_perm_kpart[bonds_iperm_kpart[bond]] = bond;

    mol->bonds.graph.renumber_vertices(bonds_perm_kpart, bonds_iperm_kpart);
  } else {
    partitioning[0] = 0;
    partitioning[1] = mol->bonds.num;
    partitioning[2] = partitioning[1];
  }

  std::vector<int> schur_solver_perm;
  schur_solver = std::unique_ptr<SchurLinearSolver>(
      new SchurLinearSolver(mol->bonds.graph, upper_tri, partitioning, schur_solver_perm));

  // Construct the final bond permutation.
  std::vector<int> bonds_perm;
  if (!bonds_perm_kpart.empty()) {
    bonds_perm.resize(mol->bonds.num);
    for (int bond = 0; bond < mol->bonds.num; ++bond)
      bonds_perm[bond] = bonds_perm_kpart[schur_solver_perm[bond]];
  } else {
    bonds_perm = std::move(schur_solver_perm);
  }
  // The bond graph is no longer needed, so do not renumber it (it has already
  // been renumbered after computing the partition).
  mol->renumber_bonds(bonds_perm, false);

  make_weights();

  // The current_lagr vectors must be as large as the rhs in order to swap them.
  current_lagr.resize(schur_solver->part_data.size());
  for (size_t p = 0; p < schur_solver->part_data.size(); ++p)
    current_lagr[p].resize(schur_solver->part_data[p].rhs.size());
}

/* ----------------------------------------------------------------------
   Compute the part of the right-hand side g(x) for rows [gstart, gend).
   Returns the largest relative (squared) bond-length violation of that part.
   The bond vectors use raw differences: the partner index already refers to
   the closest periodic image (Domain::closest_image), so no minimum-image
   correction is needed here.
------------------------------------------------------------------------- */

real Ilves::make_rhs_scalar(double **const x, double **const xprime, const bool compute_x_ab,
                            real *const rhs, const int gstart, const int gend, const int lstart)
{
  const bool xprime_ab_empty = xprime_ab.back().empty();

  real rel = 0;

  for (int grow = gstart, lrow = lstart; grow < gend; ++grow, ++lrow) {
    const int a = mol->bonds.atom1[grow];
    const int b = mol->bonds.atom2[grow];

    if (compute_x_ab)
      for (int d = 0; d < DIM; ++d) x_ab[d][grow] = x[b][d] - x[a][d];

    real rcd[DIM];
    for (int d = 0; d < DIM; ++d) rcd[d] = xprime[b][d] - xprime[a][d];
    if (!xprime_ab_empty)
      for (int d = 0; d < DIM; ++d) xprime_ab[d][grow] = rcd[d];

    const real scalar = rcd[XX] * rcd[XX] + rcd[YY] * rcd[YY] + rcd[ZZ] * rcd[ZZ];

    rhs[lrow] = 0.5 * (scalar - mol->bonds.sigma2[grow]);

    rel = std::max(rel, std::abs(rhs[lrow]) * mol->bonds.invsigma2[grow]);
  }

  return rel;
}

real Ilves::make_rhs(const int partition, double **const x, double **const xprime,
                     const bool compute_x_ab)
{
  auto &pdata = schur_solver->part_data[partition];

  // Nullify the schur entries of the rhs.
  for (int row = pdata.local_rows; row < pdata.local_rows + pdata.schur_rows; ++row)
    pdata.rhs[row] = 0;

  return make_rhs_scalar(x, xprime, compute_x_ab, pdata.rhs.data(), pdata.part[0], pdata.part[1], 0);
}

/* ----------------------------------------------------------------------
   Compute the left-hand side (Jacobian) of partition PARTITION from the bond
   vectors xab1 and xab2 (each x_ab or xprime_ab).
------------------------------------------------------------------------- */

void Ilves::make_lhs_scalar(const int partition, const BondVecs &xab1, const BondVecs &xab2,
                            const int lrowstart)
{
  auto &pdata = schur_solver->part_data[partition];
  const auto &weights = part_lhs_weights[partition];

  const int lrowend = pdata.lhs.size();
  for (int lrow = lrowstart; lrow < lrowend; ++lrow) {
    const int grow = pdata.grows[lrow];
    const int gcol = pdata.gcols[lrow];

    const real scalar = xab1[XX][grow] * xab2[XX][gcol] + xab1[YY][grow] * xab2[YY][gcol] +
        xab1[ZZ][grow] * xab2[ZZ][gcol];

    pdata.lhs[lrow] = weights[lrow] * scalar;
  }
}

void Ilves::make_lhs(const int partition, const BondVecs &xab1, const BondVecs &xab2)
{
  make_lhs_scalar(partition, xab1, xab2, 0);
}

/* ---------------------------------------------------------------------- */

void Ilves::update_current_lagr(const int partition, const bool first_time)
{
  auto &pdata = schur_solver->part_data[partition];

  if (first_time) {
    std::swap(current_lagr[partition], pdata.rhs);
  } else {
    auto *rhs_data = pdata.rhs.data();
    auto *current_lagr_data = current_lagr[partition].data();
    for (int lrow = 0; lrow < pdata.local_rows; ++lrow) current_lagr_data[lrow] += rhs_data[lrow];
  }
}

/* ----------------------------------------------------------------------
   Accumulate this Newton iteration's position increments into dx, for both
   atoms of every constraint (the atoms may be local or ghost).  dx must be
   pre-zeroed; the caller reverse-sums dx to the owning ranks and applies it.
------------------------------------------------------------------------- */

void Ilves::add_global_virial(double *const v6, const real inv_dtfsq) const
{
  // for constraint k the force on atom a is +lambda_k*inv_dtfsq*r_k (r_k = x_b -
  // x_a) and on atom b is the negative of that; the pairwise virial contribution
  // is (x_a - x_b) (x) f_a = -lambda_k*inv_dtfsq * r_k (x) r_k.
  for (size_t p = 0; p < schur_solver->part_data.size(); ++p) {
    const auto &pdata = schur_solver->part_data[p];
    for (int grow = pdata.part[0], lrow = 0; grow < pdata.part[1]; ++grow, ++lrow) {
      const real s = -current_lagr[p][lrow] * inv_dtfsq;
      const real rx = x_ab[XX][grow], ry = x_ab[YY][grow], rz = x_ab[ZZ][grow];
      v6[0] += s * rx * rx;
      v6[1] += s * ry * ry;
      v6[2] += s * rz * rz;
      v6[3] += s * rx * ry;
      v6[4] += s * rx * rz;
      v6[5] += s * ry * rz;
    }
  }
}

/* ---------------------------------------------------------------------- */

real Ilves::recompute(double **const x, double **const xprime, const bool first_iter)
{
  update_current_lagr(0, first_iter);
  return make_rhs(0, x, xprime, false);
}

/* ---------------------------------------------------------------------- */

void Ilves::accumulate_increment(const int partition, double **const dx) const
{
  const auto &pdata = schur_solver->part_data[partition];

  for (int grow = pdata.part[0], lrow = 0; grow < pdata.part[1]; ++grow, ++lrow) {
    const int a = mol->bonds.atom1[grow];
    const int b = mol->bonds.atom2[grow];

    const real rhs_a = pdata.rhs[lrow] * mol->atoms.invmass[a];
    const real rhs_b = pdata.rhs[lrow] * mol->atoms.invmass[b];
    for (int d = 0; d < DIM; ++d) {
      dx[a][d] += rhs_a * x_ab[d][grow];
      dx[b][d] -= rhs_b * x_ab[d][grow];
    }
  }
}

/* ---------------------------------------------------------------------- */

bool Ilves::disjoint_mol(const int submol_max_size) const
{
  auto &graph = mol->bonds.graph;

  std::vector<bool> visited(graph.num_nodes(), false);

  for (int bond = 0; bond < graph.num_nodes(); ++bond) {
    if (visited[bond]) continue;

    visited[bond] = true;
    int nneighs = 1;

    std::function<void(const int)> count_neighs = [&](const int node) -> void {
      for (int k = graph.xadj[node]; k < graph.xadj[node + 1]; ++k) {
        const int neigh = graph.adj[k];
        if (visited[neigh]) continue;
        visited[neigh] = true;
        ++nneighs;
        if (nneighs > submol_max_size) return;
        count_neighs(neigh);
      }
    };

    count_neighs(bond);

    if (nneighs > submol_max_size) return false;
  }

  return true;
}

/* ----------------------------------------------------------------------
   Precompute the weights of the entries of the lhs.  See the extended
   derivation in the GROMACS reference (CCKM); the weights fold in the masses
   and signs so the Jacobian A = Dg(r) inv(M) Dg(s)' can be assembled one entry
   at a time.
------------------------------------------------------------------------- */

void Ilves::make_weights()
{
  const int nparts = schur_solver->part_data.size();
  part_lhs_weights.resize(nparts);

  for (int p = 0; p < nparts; ++p) {
    const auto &pdata = schur_solver->part_data[p];
    const auto &pmatrix = pdata.fill_matrix;

    auto &pweights = part_lhs_weights[p];
    pweights.resize(pmatrix.num_edges());

    for (int i = 0; i < pmatrix.num_edges(); ++i) {
      const int row = pdata.grows[i];
      const int col = pdata.gcols[i];

      const int arow1 = mol->bonds.atom1[row];
      const int arow2 = mol->bonds.atom2[row];

      if (pdata.is_fillin[i]) {
        pweights[i] = 0;
      } else if (row != col) {
        const int acol1 = mol->bonds.atom1[col];
        const int acol2 = mol->bonds.atom2[col];

        const int common = ((arow1 == acol1) || (arow1 == acol2)) ? arow1 : arow2;

        pweights[i] = mol->atoms.invmass[common];

        if ((arow1 == acol2) || (arow2 == acol1)) pweights[i] = -pweights[i];
      } else {
        pweights[i] = mol->atoms.invmass[arow1] + mol->atoms.invmass[arow2];
      }
    }
  }
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
