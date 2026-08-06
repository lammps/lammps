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
   ILVES constraint solver: sparse direct (LU) solver.
   Specialized from the GROMACS 2021 ILVES Schur-complement solver (LGPL-2.1),
   src/gromacs/mdlib/schur_linear_solver.cpp.  Serial single-block port: the
   Schur-complement partitioning reduces to a plain sparse LU solve.  The
   minimal-degree reordering / fill-in generator is retained from upstream.
   See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#include "ilves_solver.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <memory_resource>
#include <tuple>
#include <utility>
#include <vector>

#include "ilves_graph.h"

namespace LAMMPS_NS::ILVES {

SparseDirectSolver::SparseDirectSolver(Graph &matrix, std::vector<int> &perm)
{
  nrows = matrix.num_nodes();

  FillMatrixGenerator fill_matrix_generator(matrix);
  // bind plain reference variables to the returned tuple (rather than a
  // structured binding) to keep the names usable with the widest set of C++
  // compilers/standards
  auto fill_tuple = fill_matrix_generator.get_fill_matrix();
  Graph &gfill_matrix = std::get<0>(fill_tuple);
  std::vector<bool> &gis_fillin = std::get<1>(fill_tuple);
  perm = std::move(std::get<2>(fill_tuple));

  populate(gfill_matrix, gis_fillin);
}

void SparseDirectSolver::populate(const Graph &gfill_matrix, const std::vector<bool> &gis_fillin)
{
  const int total_entries = gfill_matrix.num_edges();

  grows.resize(total_entries);
  gcols.resize(total_entries);
  is_fillin.resize(total_entries);
  diag.resize(nrows);

  fill_matrix.nnodes = nrows;
  fill_matrix.xadj.resize(nrows + 1);
  fill_matrix.adj.resize(total_entries);

  lhs.resize(total_entries);
  rhs.resize(nrows);
  scratch.resize(nrows, 0);

  // Copy the global fill-in matrix into the solver's CSR structure.  With a
  // single block the local and global row/column indices coincide.
  int lentry = 0;
  fill_matrix.xadj[0] = 0;
  for (int row = 0; row < nrows; ++row) {
    for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
      const int col = gfill_matrix.adj[k];

      grows[lentry] = row;
      gcols[lentry] = col;

      if (row == col) { diag[row] = lentry; }

      is_fillin[lentry] = gis_fillin[k];
      fill_matrix.adj[lentry] = col;

      ++lentry;
    }
    fill_matrix.xadj[row + 1] = lentry;
  }
}

void SparseDirectSolver::LU_factor()
{
  // Isolate the adjacency lists
  const auto &adj = fill_matrix.adj;
  // Isolate the list of row indices
  const auto &xadj = fill_matrix.xadj;

  // Loop over the columns of the matrix
  for (int j = 0; j < nrows; j++) {

    // Isolate the diagonal entry A(j,j)
    const double pivot = lhs[diag[j]];
    const double invpivot = 1.0 / pivot;

    // Process all *relevant* rows below row j
    for (int r = diag[j] + 1; r < xadj[j + 1]; ++r) {
      // Isolate the row index each relevant row
      const int row = adj[r];

      // Expand A(row,:) into the scratch space
      for (int s = xadj[row]; s < xadj[row + 1]; ++s) { scratch[adj[s]] = lhs[s]; }

      // Compute the multiplier
      scratch[j] = scratch[j] * invpivot;

      // Do the linear update
      for (int t = diag[j] + 1; t < xadj[j + 1]; ++t) {
        scratch[adj[t]] = scratch[adj[t]] - scratch[j] * lhs[t];
      }

      // Compress scratch back into A(row,:)
      // Also, Ensure that scratch is zero
      for (int t = xadj[row]; t < xadj[row + 1]; ++t) {
        lhs[t] = scratch[adj[t]];
        scratch[adj[t]] = 0;
      }
    }
  }
}

void SparseDirectSolver::LU_forward()
{
  const auto &adj = fill_matrix.adj;
  const auto &xadj = fill_matrix.xadj;

  // Loop over the rows, removing the already-solved contributions from the rhs
  for (int i = 0; i < nrows; i++) {
    for (int k = xadj[i]; adj[k] < i; ++k) { rhs[i] -= lhs[k] * rhs[adj[k]]; }
  }
}

void SparseDirectSolver::LU_backward()
{
  const auto &adj = fill_matrix.adj;
  const auto &xadj = fill_matrix.xadj;

  // Loop backwards over the rows
  for (int i = nrows - 1; i != -1; --i) {
    // Remove the contributions from all variables with index higher than i
    for (int k = diag[i] + 1; k < xadj[i + 1]; ++k) { rhs[i] -= lhs[k] * rhs[adj[k]]; }
    // Do the central division
    rhs[i] /= lhs[diag[i]];
  }
}

void SparseDirectSolver::LU_solve()
{
  LU_forward();
  LU_backward();
}

SparseDirectSolver::FillMatrixGenerator::FillMatrixGenerator(const Graph &matrix) :
    matrix(matrix), perm(matrix.num_nodes()), iperm(matrix.num_nodes()),
    // size the initial arena buffer for the expected number of list nodes
    // (value + the two list-node link pointers), matching the up-front
    // reservation the ported growing pool used; the arena grows if needed.
    matrix_mem_pool(matrix.num_edges() * 2 * (sizeof(MatrixEntry) + 2 * sizeof(void *))),
    active_rows_mem_pool(matrix.num_nodes() * (sizeof(int) + 2 * sizeof(void *)))
{

  fillin_matrix.nnodes = matrix.num_nodes();
  fillin_matrix.xadj.resize(matrix.num_nodes() + 1);

  init_matrices();
  init_active_rows();
  compute_fillins();
  apply_permutation();

  // Compute xadj from the per-row finalized degrees, in the new ordering.
  fillin_matrix.xadj[0] = 0;
  for (int old_row = 0; old_row < fillin_matrix.nnodes; ++old_row) {
    const int row = iperm[old_row];
    fillin_matrix.xadj[row + 1] = final_matrix[old_row].size();
  }
  // Compute the prefix sum.
  for (int row = 0; row < fillin_matrix.nnodes; ++row) {
    fillin_matrix.xadj[row + 1] += fillin_matrix.xadj[row];
  }
  fillin_matrix.adj.resize(fillin_matrix.xadj.back());
  is_fillin.resize(fillin_matrix.xadj.back());

  copy_aux_to_final();
}

std::tuple<Graph &, std::vector<bool> &, std::vector<int> &>
SparseDirectSolver::FillMatrixGenerator::get_fill_matrix()
{
  return std::make_tuple(std::ref(fillin_matrix), std::ref(is_fillin), std::ref(perm));
}

void SparseDirectSolver::FillMatrixGenerator::init_matrices()
{
  const int nnodes = matrix.num_nodes();
  init_matrix.reserve(nnodes);
  final_matrix.reserve(nnodes);

  // Copy the matrix into the per-row scratch lists.
  for (int row = 0; row < nnodes; ++row) {
    std::pmr::polymorphic_allocator<MatrixEntry> row_allocator(&matrix_mem_pool);

    init_matrix.emplace_back(row_allocator);
    final_matrix.emplace_back(row_allocator);

    for (int k = matrix.xadj[row]; k < matrix.xadj[row + 1]; ++k) {
      const int col = matrix.adj[k];

      if (col == row) {
        // The diagonal stays in the final matrix from the start.
        final_matrix.back().push_back({col, false});
      } else {
        init_matrix.back().push_back({col, false});
      }
    }
  }
}

void SparseDirectSolver::FillMatrixGenerator::init_active_rows()
{
  const int nnodes = matrix.num_nodes();
  active_rows_ptrs.reserve(nnodes);

  for (int row = 0; row < nnodes; ++row) {
    // Add the row to the active rows, bucketed by its current degree.
    const int rowdeg = init_matrix[row].size();
    auto deg_list_it = active_rows.find(rowdeg);
    if (deg_list_it == active_rows.end()) {    // Does not exist yet.
      deg_list_it = active_rows
                        .emplace(std::piecewise_construct, std::forward_as_tuple(rowdeg),
                                 std::forward_as_tuple(
                                     std::pmr::polymorphic_allocator<int>(&active_rows_mem_pool)))
                        .first;
    }
    deg_list_it->second.push_front(row);

    active_rows_ptrs.push_back(deg_list_it->second.begin());
  }
}

void SparseDirectSolver::FillMatrixGenerator::compute_fillins()
{
  int pgrow = 0;

  while (!active_rows.empty()) {
    // Find the row with minimal degree.
    const int row = active_rows.begin()->second.back();

    // Update the permutation.
    perm[pgrow] = row;
    iperm[row] = pgrow;

    for (const auto [col, col_is_fillin] : init_matrix[row]) { update_neighbors(row, col); }

    // Move the remaining ids from init_matrix to final_matrix.
    const int old_deg = init_matrix[row].size();
    final_matrix[row].splice(final_matrix[row].end(), init_matrix[row]);

    // Delete the node from the active list.
    update_active_row(row, old_deg, true);

    ++pgrow;
  }
}

void SparseDirectSolver::FillMatrixGenerator::apply_permutation()
{
  const int nnodes = matrix.num_nodes();
  // Auxiliary vector for sorting.
  std::vector<MatrixEntry> sortv(nnodes);

  for (int old_row = 0; old_row < nnodes; ++old_row) {
    int nedges = 0;
    for (auto it = final_matrix[old_row].begin(); it != final_matrix[old_row].end();) {
      const int col = iperm[it->id];

      it->id = col;    // Apply the permutation.

      // Copy into the auxiliary vector.
      sortv[nedges++] = *it;
      ++it;
    }

    // Sort the elements of the row based on the new numbering.
    std::sort(sortv.begin(), sortv.begin() + nedges, [&](const auto &a, const auto &b) {
      return a.id < b.id;
    });

    // Copy sortv back to the list.
    auto it = final_matrix[old_row].begin();
    for (int i = 0; i < nedges; ++i) {
      *it = sortv[i];
      ++it;
    }
  }
}

void SparseDirectSolver::FillMatrixGenerator::copy_aux_to_final()
{
  const int nnodes = matrix.num_nodes();
  for (int old_row = 0; old_row < nnodes; ++old_row) {
    const int row = iperm[old_row];

    int edge = fillin_matrix.xadj[row];
    for (auto [col, col_is_fillin] : final_matrix[old_row]) {
      fillin_matrix.adj[edge] = col;
      is_fillin[edge] = col_is_fillin;
      ++edge;
    }
  }
}

void SparseDirectSolver::FillMatrixGenerator::update_active_row(const int row, const int old_deg,
                                                                const bool disable)
{
  const int new_deg = init_matrix[row].size();

  if (old_deg == new_deg && !disable) { return; }

  // Iterator to avoid multiple searches.
  auto old_deg_it = active_rows.find(old_deg);
  auto &old_deg_list = old_deg_it->second;

  if (disable) {
    // Remove the node from the old list.
    old_deg_list.erase(active_rows_ptrs[row]);
  } else {
    // Move the node from the old list to the new one.
    auto new_deg_list_it = active_rows.find(new_deg);
    if (new_deg_list_it == active_rows.end()) {    // Does not exist yet.
      new_deg_list_it = active_rows
                            .emplace(std::piecewise_construct, std::forward_as_tuple(new_deg),
                                     std::forward_as_tuple(std::pmr::polymorphic_allocator<int>(
                                         &active_rows_mem_pool)))
                            .first;
    }
    new_deg_list_it->second.splice(new_deg_list_it->second.end(), old_deg_list,
                                   active_rows_ptrs[row]);
    active_rows_ptrs[row] = std::prev(new_deg_list_it->second.end());
  }

  // Remove the old degree key from the map if there are no more nodes with
  // that degree.
  if (old_deg_list.empty()) { active_rows.erase(old_deg_it); }
}

void SparseDirectSolver::FillMatrixGenerator::update_neighbors(const int row, const int col)
{
  const int col_old_deg = init_matrix[col].size();

  auto row_it = init_matrix[row].begin();
  auto col_it = init_matrix[col].begin();

  while (row_it != init_matrix[row].end() && col_it != init_matrix[col].end()) {
    // Remove row from the neighbors of col.
    if (col_it->id == row) {
      // Move from init_matrix to final_matrix.
      auto col_it_tmp = col_it;
      ++col_it;
      final_matrix[col].splice(final_matrix[col].end(), init_matrix[col], col_it_tmp);
    }
    // Do not take into account the col id in the merging.
    else if (row_it->id == col) {
      ++row_it;
    }
    // Already in both lists.
    else if (row_it->id == col_it->id) {
      ++col_it;
      ++row_it;
    }
    // Already a neighbor of the column.
    else if (col_it->id < row_it->id) {
      ++col_it;
    }
    // New neighbor of the column.
    else {
      // Every new edge is a fill-in.
      init_matrix[col].insert(col_it, {row_it->id, true});
      ++row_it;
    }
  }
  // Process the remaining neighbors of col.
  // We just need to remove row from the neighbors of col.
  while (col_it != init_matrix[col].end() && col_it->id <= row) {
    if (col_it->id == row) {
      // Move from init_matrix to final_matrix.
      auto col_it_tmp = col_it;
      ++col_it;
      final_matrix[col].splice(final_matrix[col].end(), init_matrix[col], col_it_tmp);
    } else {
      ++col_it;
    }
  }
  // Process the remaining neighbors of row.
  // We append all the remaining neighbors of row to the list of neighbors
  // of col.
  while (row_it != init_matrix[row].end()) {
    // Do not take into account the neighbor id.
    if (row_it->id != col) {
      // Every new edge is a fill-in.
      init_matrix[col].push_back({row_it->id, true});
    }
    ++row_it;
  }

  update_active_row(col, col_old_deg, false);
}

double SparseDirectSolver::memory_usage() const
{
  double bytes = fill_matrix.memory_usage();
  bytes += (double) (grows.size() + gcols.size() + diag.size()) * sizeof(int);
  bytes += (double) is_fillin.size() / 8.0;    // std::vector<bool> is bit-packed
  bytes += (double) (lhs.size() + rhs.size() + scratch.size()) * sizeof(double);
  return bytes;
}
}    // namespace LAMMPS_NS::ILVES
