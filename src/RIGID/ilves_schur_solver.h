/* -*- c++ -*- ----------------------------------------------------------
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
   src/gromacs/mdlib/schur_linear_solver.{cpp,h}.  This port is serial, so the
   matrix is a single block (no MPI/OpenMP partitioning) and the Schur-complement
   layer reduces to a plain sparse LU factorization with forward/backward solve.
   The symbolic phase (minimal-degree reordering and fill-in computation) is
   retained from upstream.  See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_SCHUR_SOLVER_H
#define LMP_ILVES_SCHUR_SOLVER_H

#include <list>
#include <map>
#include <tuple>
#include <vector>

#include "ilves_graph.h"
#include "ilves_mempool.h"

/**
 * A sparse direct solver for a structurally symmetric matrix.  The matrix is
 * stored in CSR form with a fill-reducing reordering computed once at
 * construction; it is factored in place (LU) and solved by forward/backward
 * substitution.
 */

namespace LAMMPS_NS {
namespace ILVES {

class SchurLinearSolver {
public:
    // reordered adjacency (CSR) including fill-in.  the numeric factor is held
    // in lhs, indexed in lockstep with fill_matrix.adj
    Graph fill_matrix;

    // global row and column index of each stored entry (in CSR order)
    std::vector<int> grows;
    std::vector<int> gcols;

    // is entry i a fill-in (true) or an original nonzero (false)?
    std::vector<bool> is_fillin;

    // index into lhs / fill_matrix.adj of the diagonal entry of each row
    std::vector<int> diag;

    // numeric left-hand side (one per entry) and right-hand side (one per row)
    std::vector<double> lhs;
    std::vector<double> rhs;

    // scratch row used during the factorization
    std::vector<double> scratch;

    // number of rows (= number of constraints)
    int nrows;

    /**
     * Constructs the solver for a structurally symmetric adjacency MATRIX.
     * Computes a fill-reducing reordering and the resulting fill-in structure,
     * then sizes the factor / solve work arrays.
     *
     * @param matrix Structurally symmetric adjacency matrix.
     * @param perm The array will be overwritten with the permutation applied
     * to the original matrix to reduce the number of fillins. The permutation
     * is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     */
    SchurLinearSolver(Graph &matrix, std::vector<int> &perm);

    /**
     * Performs the in-place LU factorization of the matrix.  The lhs must be
     * filled before calling this function.
     */
    void LU_factor();

    /**
     * Performs the forward + backward substitution.  LU_factor must have been
     * called and the rhs filled before calling this function.
     */
    void LU_solve();

    /**
     * Estimate the memory used by the solver.
     *
     * @return The size of the solver storage in bytes.
     */
    double memory_usage() const;

private:
    /**
     * Populate the CSR and work arrays from the global fill-in matrix.
     *
     * @param gfill_matrix The global fill-in adjacency matrix.
     * @param gis_fillin For each edge in gfill_matrix, is it a fill-in?
     */
    void populate(const Graph &gfill_matrix, const std::vector<bool> &gis_fillin);

    /** Forward substitution (lower-triangular solve), in place on rhs. */
    void LU_forward();

    /** Backward substitution (upper-triangular solve), in place on rhs. */
    void LU_backward();

    // Auxiliary class to generate the minimum fill-in matrix.
    class FillMatrixGenerator {
    public:
        /**
         * Computes the minimal degree reordering of MATRIX. Also computes the
         * fillin matrix of the reordered matrix.
         *
         * @param matrix Adjacency matrix of a structurally symmetric matrix.
         * @param parts A vector of size N + 1 that contains the partitioning of
         * the matrix. p[0] is the first row of partition 0. p[1] is the first
         * row of partition 1, and the last row + 1 of partition 0.
         */
        FillMatrixGenerator(const Graph &matrix,
                            const std::vector<int> &parts);

        /**
         * Returns the fillin matrix computed in the constructor.
         *
         * @return A tuple with three elements:
         * 1. Reference to the fillin matrix.
         * 2. Reference to a vector that describes if each edge of the fillin
         * matrix, is a fillin (true) or not (false).
         * 3. Reference to the permutation applied to the original matrix to
         * reduce the number of fillins. The permutation is given as in MATLAB.
         * Example: p = [2, 1, 0] Means that Old position 2 is now position 0
         * Old position 1 is now position 1
         * Old position 0 is now position 2
         */
        std::tuple<Graph &, std::vector<bool> &, std::vector<int> &> get_fill_matrix();

    private:
        const Graph &matrix;
        const std::vector<int> &parts;

        // The fillin matrix.
        Graph fillin_matrix;
        // The new ordering of the rows.
        std::vector<int> perm;
        // The inverse permutation
        std::vector<int> iperm;
        // True if edge i is a fillin.
        std::vector<bool> is_fillin;

        // One PartitionData instance per partition.
        class PartitionData {
        public:
            FillMatrixGenerator &parent;
            const int part;   // Partition id.

            struct MatrixEntry {
                int id;
                bool is_fillin;
            };

            // Memory pools.
            std::vector<GrowingMemPool> matrix_mem_pools;
            GrowingMemPool active_rows_mem_pool;

            // As many lists as local rows.
            // These lists will hold the initial matrix and will be used as a
            // scratch for the algorithm.
            std::vector<std::list<MatrixEntry, GrowingAllocator<MatrixEntry>>> init_matrix;
            // As many lists as local rows.
            // These lists will hold the final fillin matrix previous
            // copying it to the global fillin matrix strcuture.
            std::vector<std::list<MatrixEntry, GrowingAllocator<MatrixEntry>>> final_matrix;

            // Active rows for each partition.
            // Key = degree of the rows.
            // Value = List of rows with Key degree.
            std::map<int, std::list<int, GrowingAllocator<int>>> active_rows;
            // Iterator to where the row element is in the active_rows list.
            // Used for fast removal.
            std::vector<std::list<int, GrowingAllocator<int>>::iterator> active_rows_ptrs;

            /**
             * Get the first row of the partition.
             *
             * @return int The first row of the partition.
             */
            int first_row() const;

            /**
             * Get the last row + 1 of the partition.
             *
             * @return The last row + 1 of the partition.
             */
            int last_row_plus1() const;

            /**
             * Get the number of rows assigned to the partition.
             *
             * @return int The number of rows assigned to the partition.
             */
            int num_prows() const;

            /**
             * Get the number of edges assigned to the partition.
             *
             * @return int The number of edges assigned to the partition.
             */
            int num_pedges() const;

            /**
             * Is this partition the shared partition (last partition)?
             *
             * @return bool True if this partition is the shared partition.
             */
            bool ami_shared_part() const;

            /**
             * Get the local row id of a global row id. Undefined behavior if
             * the global row id does not belong to the partition.
             *
             * @param grow The global row id.
             * @return int The local row id.
             */
            int grow_to_prow(int grow) const;

            /**
             * Move the row (local id) to the corresponding key of active_rows,
             * taking into account that the old degree of the row was old_deg.
             *
             * @param lrow Local row id.
             * @param old_deg Old degree (key) of the row.
             * @param disable Remove the row from active_rows, do not move from
             * the old key to the new key.
             */
            void update_active_row(int lrow, int old_deg, bool disable);

            /**
             * Does updating this row (global id) from the current local
             * partition write into the shared partition?  (In the original
             * parallel code this was also the condition for taking a lock.)
             *
             * @param row The global id of the row.
             * @return true if the update is redirected to the shared partition.
             */
            bool redirect_to_shared(int row) const;

            /**
             * Merge the neighbors of the column with the neighbors of the row,
             * removing the row from the neighbors of the column and do not take
             * into account the column id in the merging.
             *
             * @param row The global id of the row. This row should be assigned
             * to the partition.
             * @param col The global id of the col.
             */
            void update_neighbors(int row, int col);

            /**
             * Construct a PartitionData object given the parent
             * FillMatrixGenerator and the partition id.
             *
             * @param parent Reference to the parent FillMatrixGenerator.
             * @param part Partition id.
             */
            PartitionData(FillMatrixGenerator &parent, int part);

            /**
             * Initialize init_matrix and final_matrix and the corresponding
             * memory pool.
             *
             */
            void init_matrices();

            /**
             * Initialize active_rows and the corresponding memory pool.
             *
             */
            void init_active_rows();

            /**
             * Compute the fillin matrix of the partition.
             *
             */
            void compute_fillins();

            /**
             * Update the column ids of the fillin matrix with the permutation.
             *
             */
            void apply_permutation();

            /**
             * Copy the fillin matrix of the partition (final_matrix) to the
             * global fillin matrix (parent.fillin_matrix).
             *
             */
            void copy_aux_to_final();
        };

        // One PartitionData instance per partition.
        std::vector<PartitionData> part_data;
    };
};
}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif
