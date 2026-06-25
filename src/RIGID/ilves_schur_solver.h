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
   ILVES constraint solver: parallel Schur-complement sparse direct solver.
   Ported near-verbatim from GROMACS 2021 ILVES (LGPL-2.1),
   src/gromacs/mdlib/schur_linear_solver.{cpp,h}.  Only the namespace, the
   include paths, and the real/AlignedAllocator/OpenMP spellings (supplied by
   ilves_compat.h) differ from upstream.  See ilves_graph.h for attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_SCHUR_SOLVER_H
#define LMP_ILVES_SCHUR_SOLVER_H

#include <list>
#include <map>
#include <set>
#include <vector>

#include "ilves_compat.h"
#include "ilves_graph.h"
#include "ilves_mempool.h"

/**
 * A class that can be used to solve linear systems of equations in parallel
 * using the Schur complement method. This class only works with structurally
 * symmetric matrices.
 */

namespace LAMMPS_NS {
namespace ILVES {

class SchurLinearSolver {
public:
    class PartitionData {
    public:
        int part[2];      // Global indices assigned to this partition.
                          // part[0] = first index of the partition.
                          // part[1] = last index of the partition + 1.
        int local_rows;   // The number of rows in the local subsystem.
        int schur_rows;   // The number of rows in the Schur complement
                          // subsystem.

        Graph fill_matrix;

        // Local rows and colums to global (whole matrix) rows and columns.
        std::vector<int, AlignedAllocator<int>> grows;
        std::vector<int, AlignedAllocator<int>> gcols;

        // Is entry i a fillin?
        std::vector<bool> is_fillin;

        // The index of the diagonal entries.
        std::vector<int> diag;

        // Mapping from local Schur idx to shared Schur block (last block) idx.
        struct localToSharedSchurMap {
            int local_idx;
            int schur_idx;
        };

        // Left-hand-side and right-hand-side local to shared Schur maps.
        std::vector<localToSharedSchurMap> lhs_local_to_shared;
        std::vector<localToSharedSchurMap> rhs_local_to_shared;

        // Left-hand-side and right-hand-side data.
        std::vector<real, AlignedAllocator<real>> lhs;
        std::vector<real, AlignedAllocator<real>> rhs;

        // Scratch vector used during the factorization.
        std::vector<real, AlignedAllocator<real>> scratch;

        /**
         * Populates the partition data given the global fill-in matrix,
         * the global fill-in, the partition id and the partitioning.
         *
         * @param gfill_matrix The global fill-in adjacency matrix.
         * @param gis_fillin For each edge in the global fill-in matrix, is it a
         * fill-in (true) or not (false)?
         * @param parts A vector that contains the partitioning of the
         * matrix. p[0] is the first row of partition 0. p[1] is the first row
         * of partition 1, and the last row + 1 of partition 0.
         * ...
         * @param part_id The id of the partition.
         */
        void populate_part(const Graph &gfill_matrix,
                           const std::vector<bool> &gis_fillin,
                           const std::vector<int> &parts,
                           int part_id);

        /**
         * Performs the LU factorization of the local partition.
         * The factorization is performed in-place. The local left-hand-side
         * of the system must be overwritten previos calling this function.
         *
         */
        void LU_factor();

        /**
         * Performs the forward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function LU_factor must be called and the rhs of the partitions must
         * be overwritten.
         *
         */
        void LU_forward();

        /**
         * Performs the backward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function LU_factor and LU_forward must be called.
         *
         */
        void LU_backward();

        /**
         * Performs the Cholesky factorization of the local partition.
         * The factorization is performed in-place. The local left-hand-side
         * of the system must be overwritten previos calling this function.
         *
         */
        void cholesky_factor();

        /**
         * Performs the forward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function cholesky_factor must be called and the rhs of the partitions
         * must be overwritten.
         *
         */
        void cholesky_forward();

        /**
         * Performs the backward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function cholesky_factor and cholesky_forward must be called.
         *
         */
        void cholesky_backward();

        /**
         * Performs the LDLT factorization of the local partition.
         * The factorization is performed in-place. The local left-hand-side
         * of the system must be overwritten previos calling this function.
         *
         */
        void LDLT_factor();

        /**
         * Performs the forward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function LDLT_factor must be called and the rhs of the partitions
         * must be overwritten.
         *
         */
        void LDLT_forward();

        /**
         * Performs the backward substitution of the local partition.
         * The substitution is performed in-place. Prior to calling this
         * function LDLT_factor and LDLT_forward must be called.
         *
         */
        void LDLT_backward();

    private:
        /**
         * Populates the local partition data given the global fill-in matrix,
         * the global fill-in vector and the rows of the Schur partition.
         *
         * @param gfill_matrix The global fill-in adjacency matrix.
         * @param gis_fillin For each edge in the global fill-in matrix, is it a
         * fill-in (true) or not (false)?
         * @param part_schur part_schur[0] = first row of the Schur partition.
         *                   part_schur[1] = last row of the Schur partition
         * + 1.
         */
        void populate_local_part(const Graph &gfill_matrix,
                                 const std::vector<bool> &gis_fillin,
                                 const int part_schur[2]);

        /**
         * Populates the Schur partition data given the global fill-in matrix,
         * the global fill-in vector.
         *
         * @param gfill_matrix The global fill-in adjacency matrix.
         * @param gis_fillin For each edge in the global fill-in matrix, is it a
         * fill-in (true) or not (false)?
         */
        void populate_schur_part(const Graph &gfill_matrix,
                                 const std::vector<bool> &gis_fillin);
    };

    // The data of each partition. The last entry of the vector is the Schur
    // (shared) partition.
    std::vector<PartitionData> part_data;

    /**
     * Constructs a SchurLinearSolver object given an adjacency matrix and
     * n partitions. The matrix must be structurally symmetric, and the rows of
     * the first n-1 partitions only can have non-zero entries in their own
     * columns or columns of the last partition. Example of a valid partitioned
     * adjacency matrix:
     *
     *   ------------
     *   |xxx|  | x |
     *   |xxf|  | f |
     *   |xfx|  | fx|
     *   |-----------
     *   |   |x |xx |
     *   |   | x| xx|
     *   |-----------
     *   |   |x |xfx|
     *   |xff|xx|fxf|
     *   |  x| x|xfx|
     *   ------------
     *
     * Partition 0 has rows 0-2.
     * Partition 1 has rows 3-4.
     * Partition 2 has rows 5-7.
     *
     * In partition 0 there are non-zero columns between columns 3-4.
     * In partition 1 there are non-zero columns between columns 0-2.
     * Partition 2 is connected to any partition.
     *
     * The constructor uses OpenMP to speedup its execution. The performance
     * is optimal when there are N - 1 available threads.
     *
     * @param matrix Structurally symmetric adjacency matrix.
     * @param upper_tri If true, the solver will only use the upper triangular
     * part of the matrix. Set this parameter to true if you want to use
     * the Cholesky or LDLT factorizations. Set it to false if you want to use
     * the LU factorization.
     * @param parts A vector of size N + 1 that contains the partition of the
     * matrix. p[0] is the first row of partition 0. p[1] is the first row of
     * partition 1, and the last row + 1 of partition 0.
     * ...
     * @param perm The array will be overwritten with the permutation applied
     * to the original matrix to reduce the number of fillins. The permutation
     * is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     */
    SchurLinearSolver(Graph &matrix,
                      bool upper_tri,
                      const std::vector<int> &parts,
                      std::vector<int> &perm);

    /**
     * Performs the parallel LU factorization of the linear-system described by
     * the object. Prior to calling this function the lhs of the partitions must
     * be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void LU_factor();

    /**
     * Performs the parallel forward+backward substitution of linear-system
     * described by the object. Prior to calling this function LU_factor must
     * be called and the rhs of the partitions must be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void LU_solve();

    /**
     * Performs the parallel Cholesky factorization of the linear-system
     * described by the object. Prior to calling this function the lhs of the
     * partitions must be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void cholesky_factor();

    /**
     * Performs the parallel forward+backward substitution of linear-system
     * described by the object. Prior to calling this function cholesky_factor
     * must be called and the rhs of the partitions must be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void cholesky_solve();

    /**
     * Performs the parallel LDLT factorization of the linear-system
     * described by the object. Prior to calling this function the lhs of the
     * partitions must be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void LDLT_factor();

    /**
     * Performs the parallel forward+backward substitution of linear-system
     * described by the object. Prior to calling this function LDLT_factor must
     * be called and the rhs of the partitions must be overwritten.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     */
    void LDLT_solve();

private:
    /**
     * Helper function to apply a factorization to the linear system using
     * one of the available PartitionData factorization functions.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     * @param factor_function Pointer to the desired PartitionData factorization
     * function.
     */
    void factor(void (PartitionData::*factor_function)());

    /**
     * Helper function to apply the forward+backward substitution to the linear
     * system using one of the available PartitionData forward+backward
     * functions.
     *
     * The function uses OpenMP for parallelization. The performance is optimal
     * when there are as many available threads as number of partitions - 1.
     *
     * @param factor_function Pointer to the desired PartitionData forward
     * function.
     * @param backward_function Pointer to the desired PartitionData backward
     * function.
     */
    void solve(void (PartitionData::*forward_function)(),
               void (PartitionData::*backward_function)());

    // Auxiliary class to generate the minimum fill-in matrix.
    class FillMatrixGenerator {
    public:
        /**
         * Computes the minimal degree reordering of MATRIX. Also computes the
         * fillin matrix of the reordered matrix.
         *
         * @param matrix Adjacency matrix of a structurally symmetric matrix.
         * @param upper_tri Only keep the upper triangular part of the matrix
         * when computing the fillin matrix.
         * @param parts A vector of size N + 1 that contains the partitioning of
         * the matrix. p[0] is the first row of partition 0. p[1] is the first
         * row of partition 1, and the last row + 1 of partition 0.
         */
        FillMatrixGenerator(const Graph &matrix,
                            bool upper_tri,
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
        const bool upper_tri;
        const std::vector<int> &parts;

        // One mutex per row of the shared partition.
        std::vector<omp_lock_t> shared_rows_locks;

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
             * Should I use a lock to access the row (global id)?
             * A pointer to the lock if true, nullptr otherwise.
             *
             * @param row The global id of the row.
             * @return omp_lock_t* A pointer to the lock that should be used to
             * access the row. Nullptr if no lock should be used.
             */
            omp_lock_t *should_lock(int row) const;

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