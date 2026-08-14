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

#ifndef LMP_ILVES_SOLVER_H
#define LMP_ILVES_SOLVER_H

#include <list>
#include <map>
#include <memory_resource>
#include <tuple>
#include <vector>

#include "ilves_graph.h"

/**
 * A sparse direct solver for a structurally symmetric matrix.  The matrix is
 * stored in CSR form with a fill-reducing reordering computed once at
 * construction; it is factored in place (LU) and solved by forward/backward
 * substitution.
 */

namespace LAMMPS_NS::ILVES {

class SparseDirectSolver {
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
  SparseDirectSolver(Graph &matrix, std::vector<int> &perm);

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
  [[nodiscard]] double memory_usage() const;

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

  // Auxiliary class that computes a fill-reducing minimal-degree reordering of
  // a structurally symmetric matrix, together with the resulting fill-in.
  class FillMatrixGenerator {
   public:
    /**
         * Computes the minimal-degree reordering of MATRIX along with the
         * fill-in matrix of the reordered matrix.
         *
         * @param matrix Adjacency matrix of a structurally symmetric matrix.
         */
    FillMatrixGenerator(const Graph &matrix);

    /**
         * Returns the fill-in matrix computed in the constructor.
         *
         * @return A tuple with three elements:
         * 1. Reference to the fill-in matrix.
         * 2. Reference to a vector that flags each edge of the fill-in matrix as
         * a fill-in (true) or an original nonzero (false).
         * 3. Reference to the permutation applied to the original matrix to
         * reduce the number of fillins. The permutation is given as in MATLAB.
         * Example: p = [2, 1, 0] Means that Old position 2 is now position 0
         * Old position 1 is now position 1
         * Old position 0 is now position 2
         */
    std::tuple<Graph &, std::vector<bool> &, std::vector<int> &> get_fill_matrix();

   private:
    struct MatrixEntry {
      int id;
      bool is_fillin;
    };

    const Graph &matrix;

    // The fillin matrix.
    Graph fillin_matrix;
    // The new ordering of the rows.
    std::vector<int> perm;
    // The inverse permutation.
    std::vector<int> iperm;
    // True if edge i is a fillin.
    std::vector<bool> is_fillin;

    // Monotonic (append-only) arenas backing the per-row scratch lists and
    // the active-row lists: the symbolic phase only ever grows these lists
    // and frees them all at once, so a bump allocator with a no-op
    // deallocate is the right fit.  std::pmr::monotonic_buffer_resource is
    // the standard C++17 equivalent of the ported GROMACS growing pool.
    // Declared before the containers below so it outlives them.
    std::pmr::monotonic_buffer_resource matrix_mem_pool;
    std::pmr::monotonic_buffer_resource active_rows_mem_pool;

    // One list per row.  init_matrix holds the working adjacency (and is
    // used as scratch by the elimination); final_matrix accumulates the
    // finalized adjacency before it is copied to the global fillin matrix.
    std::vector<std::pmr::list<MatrixEntry>> init_matrix;
    std::vector<std::pmr::list<MatrixEntry>> final_matrix;

    // Active rows bucketed by degree (key = degree, value = list of rows),
    // with an iterator per row into its bucket for fast removal.
    std::map<int, std::pmr::list<int>> active_rows;
    std::vector<std::pmr::list<int>::iterator> active_rows_ptrs;

    /**
         * Move row ROW between degree buckets after its degree changed from
         * OLD_DEG (DISABLE removes it from the active set instead of moving it).
         */
    void update_active_row(int row, int old_deg, bool disable);

    /**
         * Merge the neighbors of ROW into the neighbors of COL (recording new
         * edges as fill-ins) and remove ROW from the neighbors of COL.
         */
    void update_neighbors(int row, int col);

    /** Initialize init_matrix / final_matrix from the input matrix. */
    void init_matrices();

    /** Initialize the degree-bucketed active-row structure. */
    void init_active_rows();

    /** Run the minimal-degree elimination, filling perm and the fill-in. */
    void compute_fillins();

    /** Renumber the finalized adjacency lists with the permutation. */
    void apply_permutation();

    /** Copy the finalized adjacency lists into fillin_matrix / is_fillin. */
    void copy_aux_to_final();
  };
};
}    // namespace LAMMPS_NS::ILVES

#endif
