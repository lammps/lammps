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

#include "ilves_schur_solver.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "ilves_graph.h"
#include "ilves_mempool.h"

namespace LAMMPS_NS {
namespace ILVES {

SchurLinearSolver::SchurLinearSolver(Graph &matrix, std::vector<int> &perm) {
    nrows = matrix.num_nodes();

    // a single block spanning all rows; the empty trailing Schur block keeps
    // the fill-in generator's partitioned interface unchanged
    const std::vector<int> parts = {0, nrows, nrows};

    FillMatrixGenerator fill_matrix_generator(matrix, parts);
    // bind plain reference variables to the returned tuple (rather than a
    // structured binding) to keep the names usable with the widest set of C++
    // compilers/standards
    auto fill_tuple = fill_matrix_generator.get_fill_matrix();
    Graph &gfill_matrix = std::get<0>(fill_tuple);
    std::vector<bool> &gis_fillin = std::get<1>(fill_tuple);
    perm = std::move(std::get<2>(fill_tuple));

    populate(gfill_matrix, gis_fillin);
}

void SchurLinearSolver::populate(const Graph &gfill_matrix,
                                 const std::vector<bool> &gis_fillin) {
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

            if (row == col) {
                diag[row] = lentry;
            }

            is_fillin[lentry] = gis_fillin[k];
            fill_matrix.adj[lentry] = col;

            ++lentry;
        }
        fill_matrix.xadj[row + 1] = lentry;
    }
}

void SchurLinearSolver::LU_factor() {
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
            for (int s = xadj[row]; s < xadj[row + 1]; ++s) {
                scratch[adj[s]] = lhs[s];
            }

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

void SchurLinearSolver::LU_forward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop over the rows, removing the already-solved contributions from the rhs
    for (int i = 0; i < nrows; i++) {
        for (int k = xadj[i]; adj[k] < i; ++k) {
            rhs[i] -= lhs[k] * rhs[adj[k]];
        }
    }
}

void SchurLinearSolver::LU_backward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop backwards over the rows
    for (int i = nrows - 1; i != -1; --i) {
        // Remove the contributions from all variables with index higher than i
        for (int k = diag[i] + 1; k < xadj[i + 1]; ++k) {
            rhs[i] -= lhs[k] * rhs[adj[k]];
        }
        // Do the central division
        rhs[i] /= lhs[diag[i]];
    }
}

void SchurLinearSolver::LU_solve() {
    LU_forward();
    LU_backward();
}

SchurLinearSolver::FillMatrixGenerator::FillMatrixGenerator(const Graph &matrix,
                                                            const std::vector<int> &parts)
    : matrix(matrix),
      parts(parts),
      perm(matrix.num_nodes()),
      iperm(matrix.num_nodes()) {

    const int nparts = parts.size() - 1;

    fillin_matrix.nnodes = matrix.num_nodes();
    fillin_matrix.xadj.resize(matrix.num_nodes() + 1);

    part_data.reserve(nparts);
    for (int part = 0; part < nparts; ++part) {
        part_data.emplace_back(*this, part);
    }

    {
        // Initialize shared partition
        { part_data[nparts - 1].init_matrices(); }


        for (int part = 0; part < nparts - 1; ++part) {
            // Initialize private partitions.
            part_data[part].init_matrices();

            // Compute the fillins fo the private partitions.
            part_data[part].init_active_rows();
            part_data[part].compute_fillins();
        }   // Implicit wait.

        {
            // Compute the fillins fo the private partitions.
            part_data[nparts - 1].init_active_rows();
            part_data[nparts - 1].compute_fillins();
        }


        for (int part = 0; part < nparts; ++part) {
            // Apply the permutation in all the partitions.
            part_data[part].apply_permutation();
        }   // Implicit wait.

        // Compute xadj.
        {
            fillin_matrix.xadj[0] = 0;

            // Get the degree of each row.
            for (auto &pd : part_data) {
                for (int old_row = pd.first_row(); old_row < pd.last_row_plus1();
                     ++old_row) {
                    const int row = iperm[old_row];
                    const int lold_row = pd.grow_to_prow(old_row);

                    fillin_matrix.xadj[row + 1] = pd.final_matrix[lold_row].size();
                }
            }

            // Compute the prefix sum.
            for (int row = 0; row < fillin_matrix.nnodes; ++row) {
                fillin_matrix.xadj[row + 1] += fillin_matrix.xadj[row];
            }

            fillin_matrix.adj.resize(fillin_matrix.xadj.back());
            is_fillin.resize(fillin_matrix.xadj.back());
        }


        // Partition data to the final structure
        {
            for (int part = 0; part < nparts; ++part) {
                part_data[part].copy_aux_to_final();
            }
        }
    }
}

std::tuple<Graph &, std::vector<bool> &, std::vector<int> &>
SchurLinearSolver::FillMatrixGenerator::get_fill_matrix() {
    return std::make_tuple(std::ref(fillin_matrix),
                           std::ref(is_fillin),
                           std::ref(perm));
}

SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData(FillMatrixGenerator &parent,
                                                                     const int part)
    : parent(parent), part(part), active_rows_mem_pool(num_prows()) {

    // One matrix-memory-pool per row.
    if (ami_shared_part()) {
        matrix_mem_pools.reserve(num_prows());
        for (int row = first_row(); row < last_row_plus1(); ++row) {
            // We expect more fillins in the shared partition.
            const int nedges = parent.matrix.xadj[row + 1] -
                               parent.matrix.xadj[row];
            matrix_mem_pools.emplace_back(nedges * 5);
        }
    }
    // Just one pool.
    else {
        matrix_mem_pools.emplace_back(num_pedges() * 2);
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::init_matrices() {
    init_matrix.reserve(num_prows());
    final_matrix.reserve(num_prows());

    // Copy matrix into the scratch matrices and initialize the active nodes.
    for (int row = first_row(), lrow = 0; row < last_row_plus1(); ++row, ++lrow) {
        GrowingAllocator<MatrixEntry> row_allocator(ami_shared_part()
                                                        ? &matrix_mem_pools[lrow]
                                                        : &matrix_mem_pools.back());

        init_matrix.emplace_back(row_allocator);
        final_matrix.emplace_back(row_allocator);

        for (int k = parent.matrix.xadj[row]; k < parent.matrix.xadj[row + 1]; ++k) {
            const int col = parent.matrix.adj[k];

            if (col == row) {
                // The rowth list of the final matrix just contains the row
                // itself.
                final_matrix.back().push_back({col, false});
            }
            else {
                init_matrix.back().push_back({col, false});
            }
        }
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::init_active_rows() {
    active_rows_ptrs.reserve(num_prows());

    for (int row = first_row(), lrow = 0; row < last_row_plus1(); ++row, ++lrow) {
        // Add the row to the active rows.
        const int rowdeg = init_matrix[lrow].size();
        auto deg_list_it = active_rows.find(rowdeg);
        if (deg_list_it == active_rows.end()) {   // Does not exist yet.
            deg_list_it = active_rows
                              .emplace(std::piecewise_construct,
                                       std::forward_as_tuple(rowdeg),
                                       std::forward_as_tuple(GrowingAllocator<int>(&active_rows_mem_pool)))
                              .first;
        }
        deg_list_it->second.push_front(row);

        active_rows_ptrs.push_back(deg_list_it->second.begin());
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::compute_fillins() {
    int pgrow = first_row();

    while (!active_rows.empty()) {
        // Find the row with minimal degree.
        const int row = active_rows.begin()->second.back();
        const int lrow = grow_to_prow(row);

        // Update the permutation
        parent.perm[pgrow] = row;
        parent.iperm[row] = pgrow;

        for (const auto [col, col_is_fillin] : init_matrix[lrow]) {
            update_neighbors(row, col);
        }

        // Move the remaining ids from init_matrix to final_matrix.
        const int old_deg = init_matrix[lrow].size();
        final_matrix[lrow].splice(final_matrix[lrow].end(), init_matrix[lrow]);

        // Delete the node from the active list
        update_active_row(lrow, old_deg, true);

        ++pgrow;
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::apply_permutation() {
    // Auxiliary vector for sorting.
    std::vector<MatrixEntry> sortv(parent.matrix.num_nodes());

    for (int old_row = first_row(); old_row < last_row_plus1(); ++old_row) {
        const int lold_row = grow_to_prow(old_row);

        int nedges = 0;
        for (auto it = final_matrix[lold_row].begin();
             it != final_matrix[lold_row].end();) {
            const int col = parent.iperm[it->id];

            it->id = col;   // Apply the permutation.

            // Copy into the auxiliary vector.
            sortv[nedges++] = *it;
            ++it;
        }

        // Sort the elements of the row based on the new numbering.
        std::sort(sortv.begin(),
                  sortv.begin() + nedges,
                  [&](const auto &a, const auto &b) {
                      return a.id < b.id;
                  });

        // Copy sortv back to the list.
        auto it = final_matrix[lold_row].begin();
        for (int i = 0; i < nedges; ++i) {
            *it = sortv[i];
            ++it;
        }
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::copy_aux_to_final() {
    for (int old_row = first_row(); old_row < last_row_plus1(); ++old_row) {
        const int row = parent.iperm[old_row];
        const int lold_row = grow_to_prow(old_row);

        int edge = parent.fillin_matrix.xadj[row];
        for (auto [col, col_is_fillin] : final_matrix[lold_row]) {
            parent.fillin_matrix.adj[edge] = col;
            // Remember that is_fillin is a vector<bool> and it is not
            // thread-safe.
            parent.is_fillin[edge] = col_is_fillin;
            ++edge;
        }
    }
}

bool SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData::ami_shared_part()
    const {
    return part == parent.parts.size() - 2;
}

int SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData::first_row() const {
    return parent.parts[part];
}
int SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData::last_row_plus1() const {
    return parent.parts[part + 1];
}

int SchurLinearSolver::FillMatrixGenerator::PartitionData::grow_to_prow(const int grow) const {
    return grow - first_row();
}

int SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData::num_prows() const {
    return last_row_plus1() - first_row();
}
int SchurLinearSolver::FillMatrixGenerator::PartitionData::PartitionData::num_pedges() const {
    return parent.matrix.xadj[last_row_plus1()] - parent.matrix.xadj[first_row()];
}

bool
SchurLinearSolver::FillMatrixGenerator::PartitionData::redirect_to_shared(const int row) const {
    const bool shared_part = ami_shared_part();
    const int first_shared_row = parent.part_data.back().first_row();

    return !shared_part && row >= first_shared_row;
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::
    update_active_row(const int lrow, const int old_deg, const bool disable) {

    const int new_deg = init_matrix[lrow].size();

    if (old_deg == new_deg && !disable) {
        return;
    }

    // Iterator to avoid multiple searches.
    auto old_deg_it = active_rows.find(old_deg);
    auto &old_deg_list = old_deg_it->second;

    if (disable) {
        // Remove the node from the old list.
        old_deg_list.erase(active_rows_ptrs[lrow]);
    }
    else {
        // Move the node from the old list to the new one.
        auto new_deg_list_it = active_rows.find(new_deg);
        if (new_deg_list_it == active_rows.end()) {   // Does not exist yet.
            new_deg_list_it = active_rows
                                  .emplace(std::piecewise_construct,
                                           std::forward_as_tuple(new_deg),
                                           std::forward_as_tuple(GrowingAllocator<
                                                                 int>(&active_rows_mem_pool)))
                                  .first;
        }
        new_deg_list_it->second.splice(new_deg_list_it->second.end(),
                                       old_deg_list,
                                       active_rows_ptrs[lrow]);
        active_rows_ptrs[lrow] = std::prev(new_deg_list_it->second.end());
    }

    // Remove the old degree key from the map if
    // there are no more nodes with that degree.
    if (old_deg_list.empty()) {
        active_rows.erase(old_deg_it);
    }
}

void SchurLinearSolver::FillMatrixGenerator::PartitionData::update_neighbors(const int row,
                                                                             const int col) {
    const bool shared = redirect_to_shared(col);

    auto &col_part_data = shared ? parent.part_data.back() : *this;

    const int lrow = grow_to_prow(row);
    const int lcol = col_part_data.grow_to_prow(col);

    const int col_old_deg = col_part_data.init_matrix[lcol].size();

    auto row_it = init_matrix[lrow].begin();
    auto col_it = col_part_data.init_matrix[lcol].begin();

    while (row_it != init_matrix[lrow].end() &&
           col_it != col_part_data.init_matrix[lcol].end()) {
        // Remove row from the neighbors of col.
        if (col_it->id == row) {
            // Move from init_matrix to final_matrix.
            auto col_it_tmp = col_it;
            ++col_it;
            col_part_data.final_matrix[lcol]
                .splice(col_part_data.final_matrix[lcol].end(),
                        col_part_data.init_matrix[lcol],
                        col_it_tmp);
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
            col_part_data.init_matrix[lcol].insert(col_it, {row_it->id, true});
            ++row_it;
        }
    }
    // Process the remaining neighbors of col.
    // We just need to remove row from the neighbors of col.
    while (col_it != col_part_data.init_matrix[lcol].end() && col_it->id <= row) {
        if (col_it->id == row) {
            // Move from init_matrix to final_matrix.
            auto col_it_tmp = col_it;
            ++col_it;
            col_part_data.final_matrix[lcol]
                .splice(col_part_data.final_matrix[lcol].end(),
                        col_part_data.init_matrix[lcol],
                        col_it_tmp);
        }
        else {
            ++col_it;
        }
    }
    // Process the remaining neighbors of row.
    // We append all the remaining neighbors of row the list of
    // neighbors of col.
    while (row_it != init_matrix[lrow].end()) {
        // Do not take into account the neighbor id.
        if (row_it->id != col) {
            // Every new edge is a fill-in.
            col_part_data.init_matrix[lcol].push_back({row_it->id, true});
        }
        ++row_it;
    }

    if (!shared) {
        update_active_row(lcol, col_old_deg, false);
    }
}

double SchurLinearSolver::memory_usage() const
{
    double bytes = fill_matrix.memory_usage();
    bytes += (double) (grows.size() + gcols.size() + diag.size()) * sizeof(int);
    bytes += (double) is_fillin.size() / 8.0;    // std::vector<bool> is bit-packed
    bytes += (double) (lhs.size() + rhs.size() + scratch.size()) * sizeof(double);
    return bytes;
}

}    // namespace ILVES
}    // namespace LAMMPS_NS
