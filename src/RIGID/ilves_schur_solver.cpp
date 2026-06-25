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
   ILVES constraint solver: parallel Schur-complement sparse direct solver.
   Ported near-verbatim from GROMACS 2021 ILVES (LGPL-2.1),
   src/gromacs/mdlib/schur_linear_solver.cpp.  See ilves_graph.h for full
   attribution.
------------------------------------------------------------------------- */

#include "ilves_schur_solver.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "ilves_compat.h"
#include "ilves_graph.h"
#include "ilves_mempool.h"

namespace LAMMPS_NS {
namespace ILVES {

void SchurLinearSolver::PartitionData::populate_part(const Graph &gfill_matrix,
                                                     const std::vector<bool> &gis_fillin,
                                                     const std::vector<int> &parts,
                                                     const int part_id) {

    part[0] = parts[part_id];
    part[1] = parts[part_id + 1];

    /*
     * Example of partitioned matrix. x regular entries and f for fillin
     * entries.
     *
     *  0   1  S
     * ------------
     * |xxx|  | x | 0
     * |xxf|  | f |
     * |xfx|  | fx|
     * |-----------
     * |   |x |xx | 1
     * |   | x| xx|
     * |-----------
     * |   |x |xfx| S
     * |xff|xx|fxf|
     * |  x| x|xfx|
     * ------------
     *
     * The first two vertical blocks (0, 1) are local blocks, while the
     * last block (S) is the Schur block.
     * The first two horizontal blocks (0, 1) are local blocks, while the
     * last block (S) is the Schur block.
     * Naming of blocks:
     * local-local = vertical local block and horizontal local block.
     * local-Schur = vertical local block and horizontal Schur block.
     * Schur-local = vertical Schur block and horizontal local block.
     * Schur-Schur = vertical Schur block and horizontal Schur block.
     *
     * A local partition (partitions from 0 to num_part - 2) have local-local,
     * local-Schur, Schur-local and Schur-Schur blocks. All the entries of the
     * Schur-Schur block are fillin entries.
     *
     * The last partition (num_part - 1) is the Schur partition and only has the
     * Schur-Schur block.
     *
     * const int num_part = part.size() - 1;
     */

    const int schur_part = parts.size() - 2;
    int part_schur[2] = {parts[schur_part], parts[schur_part + 1]};

    if (schur_part == part_id) {
        populate_schur_part(gfill_matrix, gis_fillin);
    }
    else {
        populate_local_part(gfill_matrix, gis_fillin, part_schur);
    }
}

void SchurLinearSolver::PartitionData::LU_factor() {
    // Isolate the adjacency lists
    const auto &adj = fill_matrix.adj;
    // Isolate the list of row indices
    const auto &xadj = fill_matrix.xadj;

    // Loop over the first local_rows columns of the matrix
    for (int j = 0; j < local_rows; j++) {

        // Isolate the diagonal entry A(j,j)
        const real pivot = lhs[diag[j]];
        const real invpivot = 1.0 / pivot;

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

void SchurLinearSolver::PartitionData::LU_forward() {
    const auto m = fill_matrix.num_nodes();   // The number of rows
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop over the m rows
    for (int i = 0; i < m; i++) {
        // Remove relevant entries from the rhs
        const int ncols = std::min(i, local_rows);
        for (int k = xadj[i]; adj[k] < ncols; ++k) {
            rhs[i] -= lhs[k] * rhs[adj[k]];
        }
    }
}

void SchurLinearSolver::PartitionData::LU_backward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop backwards over the first n rows
    for (int i = local_rows - 1; i != -1; --i) {
        // Remove the contributions from all variables with index higher than i
        for (int k = diag[i] + 1; k < xadj[i + 1]; ++k) {
            rhs[i] -= lhs[k] * rhs[adj[k]];
        }
        // Do the central division
        rhs[i] /= lhs[diag[i]];
    }
}

void SchurLinearSolver::PartitionData::cholesky_factor() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop over the first local_rows rows of the matrix
    for (int row = 0; row < local_rows; ++row) {
        // Take the square root of the diagonal entry A(row,row)
        const real pivot = std::sqrt(lhs[diag[row]]);
        const real invpivot = 1.0 / pivot;
        lhs[diag[row]] = pivot;

        /* In the next loop we scale the upperdiagonal entries of the rowth row
           with A(row,row)
              A(row,row+1:n) /= A(row,row).
           Moreover, we expand the compressed representation of the updated
           entries into the array WORK
        */
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            // Scale the subdiagonal entry
            lhs[s] *= invpivot;
            // Copy it into the auxillary vector work
            scratch[adj[s]] = lhs[s];
        }

        /* At this point we have to perform a rank 1 update of the submatrix in
           the current upper right corner, specifically

           A(row+1:n,row+1:n) -= A(row+1:n,row) * A(row+1:n,j)'

           Only the lower triangular part matters. In MATLAB we would write

             for s=row+1:n
                alpha=A(row,s);
                A(s,s:n) -= alpha*A(row,s:n);
             end

          which stresses the fact that only nonzero entries alpha = A(row,s)
          matter.

        */

        // Loop over the nonzero upperdiagonal elements of the colth row.
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            // Isolate the nonzero entry
            const real alpha = lhs[s];
            // Isolate the index of the row which we have to update
            const int col = adj[s];
            // Update the upperdiagonal part of this row
            for (int t = diag[col]; t < xadj[col + 1]; ++t) {
                lhs[t] -= alpha * scratch[adj[t]];
            }
        }

        /* At this point the jth rank 1 update is complete. We must clear the
           nonzero entries of SCRATCH, so that it does not contain garbage
           during the next iteration.

           Notice that we do not bother to fill the entire array with zeros.
           We only kill those which could be nonzero.
        */
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            scratch[adj[s]] = 0;
        }
    }
}

void SchurLinearSolver::PartitionData::cholesky_forward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop over the local_rows rows
    for (int row = 0; row < local_rows; ++row) {
        // Divide with diagonal entry
        // TODO: Tril the graph to use xadj[row] instead of diag[row].
        rhs[row] /= lhs[diag[row]];
        /* At this point x[row] has been computed and we must eliminate it from
           equations row+1, row+2, ..., n-1.
        */

        // Remove relevant entries from the rhs
        // TODO: Tril the graph to use xadj[row] instead of diag[row].
        for (int k = diag[row] + 1; k < xadj[row + 1]; ++k) {
            // Removing the influence of x[row] from row rows[diag] of the RHS
            rhs[adj[k]] -= lhs[k] * rhs[row];
        }
    }
}

void SchurLinearSolver::PartitionData::cholesky_backward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop backwards over the first local_rows rows
    for (int row = local_rows - 1; row != -1; --row) {
        // Remove the contribution of the variables x[row+1], ..., x[n-1]
        for (int k = diag[row] + 1; k < xadj[row + 1]; ++k) {
            rhs[row] -= lhs[k] * rhs[adj[k]];
        }
        // Do the central division
        rhs[row] /= lhs[diag[row]];
    }
}

void SchurLinearSolver::PartitionData::LDLT_factor() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    for (int row = 0; row < local_rows; ++row) {
        const real pivot = lhs[diag[row]];
        const real invpivot = 1.0 / pivot;

        /* In the next loop we scale the upperdiagonal entries of the rowth row
           with A(row,row)
              A(row,row+1:n) /= A(row,row).
           Moreover, we expand the compressed representation of the updated
           entries into the array SCRATCH
        */
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            // Copy it into the auxillary vector work
            scratch[adj[s]] = lhs[s];
            // Scale the subdiagonal entry AFTER copying it.
            lhs[s] *= invpivot;
        }

        // We perform an LU factorization but only using and updating the upper
        // triangular part of the matrix. Since A is symmetric, when we need to
        // read A[row][col] we can read A[col][row] instead.

        // Loop over the nonzero upperdiagonal elements of the colth row.
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            // Isolate the nonzero entry
            const real alpha = lhs[s];
            // Isolate the index of the row which we have to update
            const int col = adj[s];
            // Update the upperdiagonal part of this row
            for (int t = diag[col]; t < xadj[col + 1]; ++t) {
                lhs[t] -= alpha * scratch[adj[t]];
            }
        }

        /* At this point the jth rank 1 update is complete. We must clear the
           nonzero entries of SCRATCH, so that it does not contain garbage
           during the next iteration.

           Notice that we do not bother to fill the entire array with zeros.
           We only kill those which could be nonzero.
        */
        for (int s = diag[row] + 1; s < xadj[row + 1]; ++s) {
            scratch[adj[s]] = 0;
        }
    }
}

void SchurLinearSolver::PartitionData::LDLT_forward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop over the local_rows rows
    for (int row = 0; row < local_rows; ++row) {
        // Loop over the strictly upperdiagonal entries of the rowth row
        for (int k = diag[row] + 1; k < xadj[row + 1]; ++k) {
            rhs[adj[k]] -= lhs[k] * rhs[row];
        }
    }
}

void SchurLinearSolver::PartitionData::LDLT_backward() {
    const auto &adj = fill_matrix.adj;
    const auto &xadj = fill_matrix.xadj;

    // Loop backwards over the first n rows
    for (int row = local_rows - 1; row != -1; --row) {
        // Isolate the rowth diagonal element
        const real d = lhs[diag[row]];

        // We only are interested in the strictly upperdiagonal part of the
        // matrix.
        for (int k = diag[row] + 1; k < xadj[row + 1]; k++) {
            // Removing the influence of x[row] from row rows[diag] of the RHS
            rhs[row] -= lhs[k] * d * rhs[adj[k]];
        }

        rhs[row] /= d;
    }
}

void SchurLinearSolver::PartitionData::populate_local_part(const Graph &gfill_matrix,
                                                           const std::vector<bool> &gis_fillin,
                                                           const int part_schur[2]) {

    // Count the number of entries and compute the mapping from global cols to
    // local cols.
    std::vector<int> gcol_to_lcol(gfill_matrix.num_nodes(), -1);

    // Horizontal local blocks.
    int local_entries = 0;
    local_rows = part[1] - part[0];
    for (int row = part[0]; row < part[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.

            gcol_to_lcol[col] = 0;
            ++local_entries;
        }
    }
    // Horizontal Schur blocks.
    int schur_entries = 0;
    // Horizontal Schur-Schur blocks.
    int schur_schur_entries = 0;
    schur_rows = 0;
    for (int row = part_schur[0]; row < part_schur[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.

            /*
             * Check if entry is part of the ith Schur-local block
             * or part of the Schur-Schur block. We only want to take
             * into account rows in the schur-schur block if the rowth
             * row of the ith Schur-local is not empty.
             *
             * The marked row in the following example is not taken into
             * account in the 0th partition, as the first row of the
             * 0th local-Schur block is empty.
             *  0   1  S
             * ------------
             * |xxx|  | x | 0
             * |xxf|  | f |
             * |xfx|  | fx|
             * |-----------
             * |   |x |xx | 1
             * |   | x| xx|
             * |-----------
             * |   |x |xfx| S <--- Not taken into account in 0th partition.
             * |xff|xx|fxf|
             * |  x| x|xfx|
             * ------------
             *
             *
             */

            const bool schur_schur = part_schur[0] <= col && col < part_schur[1];
            const bool part_schur_local = part[0] <= col && col < part[1];

            // Take Schur-Schur entry into account only if there is at least one
            // entry in same column in the local-Schur block and there is at
            // least one entry in the same row in the Schur-local block.
            const bool part_schur_schur = schur_schur && gcol_to_lcol[col] != -1 &&
                                          gcol_to_lcol[row] != -1;

            if (part_schur_local || part_schur_schur) {
                gcol_to_lcol[col] = 0;
                ++schur_entries;

                if (part_schur_schur) {
                    ++schur_schur_entries;
                }
            }
        }
        // Only if current row's Schur-local block is not empty.
        if (gcol_to_lcol[row] != -1) {
            ++schur_rows;
        }
    }

    const int total_rows = local_rows + schur_rows;
    const int total_entries = local_entries + schur_entries;

    grows.resize(total_entries);
    gcols.resize(total_entries);

    is_fillin.resize(total_entries);

    diag.resize(total_rows);

    lhs_local_to_shared.resize(schur_schur_entries);
    rhs_local_to_shared.resize(schur_rows);

    fill_matrix.nnodes = total_rows;
    fill_matrix.xadj.resize(total_rows + 1);
    fill_matrix.adj.resize(total_entries);

    lhs.resize(total_entries);
    rhs.resize(total_rows);
    scratch.resize(total_rows, 0);

    // Compute the mapping from global cols to local cols.
    for (int gcol = 0, lcol = 0; gcol < gcol_to_lcol.size(); ++gcol) {
        if (gcol_to_lcol[gcol] != -1) {
            gcol_to_lcol[gcol] = lcol;
            ++lcol;
        }
    }

    // Populate the structures.

    int lentry = 0;   // Entry index.
    int lrow = 0;     // Local row index.

    fill_matrix.xadj[0] = 0;

    // Loop over the rows of the local horizontal block.
    for (int row = part[0]; row < part[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.

            // Check if entry is part of the local-local block or part of
            // the local-Schur block. No need to do this, we are always
            // in one of the two blocks.
            // const bool ith_local_local = part[part_id] <= col &&
            //                              col < part[part_id + 1];
            // const bool ith_local_schur = part[schur_part_id] <= col &&
            //                              col < part[schur_part_id + 1];

            // if (ith_local_local || ith_local_schur) {
            const int lcol = gcol_to_lcol[col];

            grows[lentry] = row;
            gcols[lentry] = col;

            if (lrow == lcol) {
                diag[lrow] = lentry;
            }

            is_fillin[lentry] = gis_fillin[k];

            fill_matrix.adj[lentry] = lcol;

            ++lentry;
        }
        fill_matrix.xadj[lrow + 1] = lentry;
        ++lrow;
    }

    // The schur_schur block local index. Used to compute the mapping.
    int schur_schur_lentry = 0;
    int lhs_local_to_shared_idx = 0;
    int rhs_local_to_shared_idx = 0;
    // Loop over the rows of the horizontal Schur block.
    for (int row = part_schur[0]; row < part_schur[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.
            const bool schur_schur = part_schur[0] <= col && col < part_schur[1];
            const bool part_schur_local = part[0] <= col && col < part[1];
            // Same logic as before.
            const bool part_schur_schur = schur_schur && gcol_to_lcol[col] != -1 &&
                                          gcol_to_lcol[row] != -1;

            if (part_schur_local || part_schur_schur) {
                const int lcol = gcol_to_lcol[col];

                grows[lentry] = row;
                gcols[lentry] = col;

                if (part_schur_local) {
                    is_fillin[lentry] = gis_fillin[k];
                }
                else if (part_schur_schur) {
                    // We need to update the mapping.
                    lhs_local_to_shared[lhs_local_to_shared_idx] = {lentry,
                                                                    schur_schur_lentry};
                    ++lhs_local_to_shared_idx;

                    // The Schur-Schur block in a local partition is always
                    // filled with fillins.
                    is_fillin[lentry] = true;
                }

                if (lrow == lcol) {
                    diag[lrow] = lentry;
                }

                fill_matrix.adj[lentry] = lcol;

                ++lentry;
            }

            if (schur_schur) {
                ++schur_schur_lentry;
            }
        }
        // Only if current row's Schur-local block is not empty.
        if (gcol_to_lcol[row] != -1) {
            const int schur_schur_row = row - part_schur[0];
            rhs_local_to_shared[rhs_local_to_shared_idx] = {lrow, schur_schur_row};
            ++rhs_local_to_shared_idx;

            fill_matrix.xadj[lrow + 1] = lentry;

            ++lrow;
        }
    }
}

void SchurLinearSolver::PartitionData::populate_schur_part(const Graph &gfill_matrix,
                                                           const std::vector<bool>
                                                               &gis_fillin) {

    local_rows = part[1] - part[0];
    schur_rows = 0;   // No schur rows in the Schur partition.

    // Count the number of entries

    int total_entries = 0;
    for (int row = part[0]; row < part[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.

            // Check if entry is part of the Schur-Schur block.
            const bool schur_schur = part[0] <= col && col < part[1];

            if (schur_schur) {
                ++total_entries;
            }
        }
    }

    // Reserve the required space.
    grows.resize(total_entries);
    gcols.resize(total_entries);

    is_fillin.resize(total_entries);

    diag.resize(local_rows);

    fill_matrix.nnodes = local_rows;
    fill_matrix.xadj.resize(local_rows + 1);
    fill_matrix.adj.resize(total_entries);

    lhs.resize(total_entries);
    rhs.resize(local_rows);
    scratch.resize(local_rows, 0);

    // Populate the structures.

    int lentry = 0;   // Entry index.
    int lrow = 0;     // Local row index.

    fill_matrix.xadj[0] = 0;

    for (int row = part[0]; row < part[1]; ++row) {
        for (int k = gfill_matrix.xadj[row]; k < gfill_matrix.xadj[row + 1]; ++k) {
            const int col = gfill_matrix.adj[k];   // Global column index.

            // Check if entry is part of the Schur-Schur block.
            const bool schur_schur = part[0] <= col && col < part[1];

            if (schur_schur) {
                // The local column index is the global column index
                // minus the first row of the Schur block.
                const int lcol = col - part[0];

                grows[lentry] = row;
                gcols[lentry] = col;

                if (lrow == lcol) {
                    diag[lrow] = lentry;
                }

                is_fillin[lentry] = gis_fillin[k];

                fill_matrix.adj[lentry] = lcol;

                ++lentry;
            }
        }
        fill_matrix.xadj[lrow + 1] = lentry;
        ++lrow;
    }
}

SchurLinearSolver::SchurLinearSolver(Graph &matrix,
                                     const bool upper_tri,
                                     const std::vector<int> &parts,
                                     std::vector<int> &perm) {

    FillMatrixGenerator fill_matrix_generator(matrix, upper_tri, parts);
    auto [fill_matrix,
          is_fillin,
          perm_aux] = std::move(fill_matrix_generator.get_fill_matrix());
    perm = std::move(perm_aux);

    const int nparts = parts.size() - 1;

    part_data.resize(nparts);

    #pragma omp parallel num_threads(nparts - 1)
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < nparts - 1; ++i) {
            part_data[i].populate_part(fill_matrix, is_fillin, parts, i);
        }

        #pragma omp master
        {
            part_data[nparts - 1].populate_part(fill_matrix,
                                                is_fillin,
                                                parts,
                                                nparts - 1);
        }
    }
}

void SchurLinearSolver::LU_factor() {
    factor(&SchurLinearSolver::PartitionData::LU_factor);
}

void SchurLinearSolver::LU_solve() {
    solve(&PartitionData::LU_forward, &PartitionData::LU_backward);
}

void SchurLinearSolver::cholesky_factor() {
    factor(&PartitionData::cholesky_factor);
}

void SchurLinearSolver::cholesky_solve() {
    solve(&PartitionData::cholesky_forward, &PartitionData::cholesky_backward);
}

void SchurLinearSolver::LDLT_factor() { factor(&PartitionData::LDLT_factor); }

void SchurLinearSolver::LDLT_solve() {
    solve(&PartitionData::LDLT_forward, &PartitionData::LDLT_backward);
}

void SchurLinearSolver::factor(void (PartitionData::*factor_function)()) {

    const bool empty_schur = part_data.back().local_rows == 0;

    /*
     * General block factorization.
     */
    const int nlocal_parts = part_data.size() - 1;
    #pragma omp for schedule(static) nowait
    for (int p = 0; p < nlocal_parts; ++p) {
        (part_data[p].*factor_function)();

        if (!empty_schur) {
            // Add the local contributions to the Schur complement LHS
            for (const auto &map : part_data[p].lhs_local_to_shared) {
                // Add the contribution
                #pragma omp atomic
                part_data.back().lhs[map.schur_idx] += part_data[p].lhs[map.local_idx];
            }
        }
    }

    if (!empty_schur) {
        /*
         * Schur block factorization.
         */
        #pragma omp barrier

        #pragma omp master
        { (part_data.back().*factor_function)(); }
    }
}

void SchurLinearSolver::solve(void (PartitionData::*forward_function)(),
                              void (PartitionData::*backward_function)()) {

    const bool empty_schur = part_data.back().local_rows == 0;

    const int nlocal_parts = part_data.size() - 1;

    /*
     * General block forward substition.
     */
    #pragma omp for schedule(static) nowait
    for (int p = 0; p < nlocal_parts; ++p) {
        (part_data[p].*forward_function)();

        if (!empty_schur) {
            for (const auto &map : part_data[p].rhs_local_to_shared) {
                #pragma omp atomic
                part_data.back().rhs[map.schur_idx] += part_data[p].rhs[map.local_idx];
            }
        }
    }

    /*
     * Schur block forward and backward substitution.
     */
    if (!empty_schur) {
        #pragma omp barrier

        #pragma omp master
        {
            // Forward sweep for the Schur complement system
            (part_data.back().*forward_function)();
            // Backward sweep for the Schur complement system
            (part_data.back().*backward_function)();
        }

        #pragma omp barrier
    }

    /*
     * General block backward substitution.
     */
    #pragma omp for schedule(static) nowait
    for (int p = 0; p < nlocal_parts; ++p) {
        if (!empty_schur) {
            // Copy back the rhs of the Schur complement.
            for (const auto &map : part_data[p].rhs_local_to_shared) {
                part_data[p].rhs[map.local_idx] = part_data.back().rhs[map.schur_idx];
            }
        }

        (part_data[p].*backward_function)();
    }
}

SchurLinearSolver::FillMatrixGenerator::FillMatrixGenerator(const Graph &matrix,
                                                            const bool upper_tri,
                                                            const std::vector<int> &parts)
    : matrix(matrix),
      upper_tri(upper_tri),
      parts(parts),
      perm(matrix.num_nodes()),
      iperm(matrix.num_nodes()) {

    const int nparts = parts.size() - 1;

    shared_rows_locks.resize(parts[nparts] - parts[nparts - 1]);
    for (auto &lock : shared_rows_locks) {
        omp_init_lock(&lock);
    }

    fillin_matrix.nnodes = matrix.num_nodes();
    fillin_matrix.xadj.resize(matrix.num_nodes() + 1);

    part_data.reserve(nparts);
    for (int part = 0; part < nparts; ++part) {
        part_data.emplace_back(*this, part);
    }

    #pragma omp parallel num_threads(nparts - 1)
    {
        // Initialize shared partition
        #pragma omp master
        { part_data[nparts - 1].init_matrices(); }

        #pragma omp barrier

        #pragma omp for schedule(static)
        for (int part = 0; part < nparts - 1; ++part) {
            // Initialize private partitions.
            part_data[part].init_matrices();

            // Compute the fillins fo the private partitions.
            part_data[part].init_active_rows();
            part_data[part].compute_fillins();
        }   // Implicit wait.

        #pragma omp master
        {
            // Compute the fillins fo the private partitions.
            part_data[nparts - 1].init_active_rows();
            part_data[nparts - 1].compute_fillins();
        }

        #pragma omp barrier

        #pragma omp for schedule(static)
        for (int part = 0; part < nparts; ++part) {
            // Apply the permutation in all the partitions.
            part_data[part].apply_permutation();
        }   // Implicit wait.

        // Compute xadj.
        #pragma omp master
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

        #pragma omp barrier

        // std::vector<bool> is not thread-safe for write accesses, since it
        // uses a bitmap. We serialize this loop as it is not performance
        // critical.
        // #pragma omp for schedule(static)
        // Partition data to the final structure
        #pragma omp master
        {
            for (int part = 0; part < nparts; ++part) {
                part_data[part].copy_aux_to_final();
            }
        }
    }

    for (auto &lock : shared_rows_locks) {
        omp_destroy_lock(&lock);
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
        const int row = parent.iperm[old_row];
        const int lold_row = grow_to_prow(old_row);

        int nedges = 0;
        for (auto it = final_matrix[lold_row].begin();
             it != final_matrix[lold_row].end();) {
            const int col = parent.iperm[it->id];

            it->id = col;   // Apply the permutation.

            // Remove element.
            if (parent.upper_tri && col < row) {
                it = final_matrix[lold_row].erase(it);
            }
            else {
                // Copy into the auxiliary vector.
                sortv[nedges++] = *it;
                ++it;
            }
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

omp_lock_t *
SchurLinearSolver::FillMatrixGenerator::PartitionData::should_lock(const int row) const {
    const bool shared_part = ami_shared_part();
    const int first_shared_row = parent.part_data.back().first_row();

    if (!shared_part && row >= first_shared_row) {
        return &parent.shared_rows_locks[row - first_shared_row];
    }
    else {
        return nullptr;
    }
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
    auto *lock = should_lock(col);

    if (lock) {
        omp_set_lock(lock);
    }

    auto &col_part_data = lock ? parent.part_data.back() : *this;

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

    if (lock) {
        omp_unset_lock(lock);
    }
    else {
        update_active_row(lcol, col_old_deg, false);
    }
}
}    // namespace ILVES
}    // namespace LAMMPS_NS
