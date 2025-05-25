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
   Contributing authors: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */


#include <iostream>
#include <iomanip>
#include "Kokkos_Core.hpp"
#include "KokkosBlas1_dot.hpp"
#include "KokkosBatched_CG.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

static constexpr double QSUMSMALL = 0.00001;

//#include <iostream>
//#include <iomanip>

template<typename ViewType>
void print_1d_view(const ViewType& view, const char* label, int max_elements = 10) {
    auto h_view = Kokkos::create_mirror_view(view);
    Kokkos::deep_copy(h_view, view);
    
    //printf("%s[%d]: ", label, (int)view.extent(0));
    printf("%s: ", label);
    int n = std::min(max_elements, (int)view.extent(0));
    for (int i = 0; i < n; i++) {
        printf("%8.6f ", (double)h_view(i));
    }
    //if (view.extent(0) > max_elements) printf("...");
    printf("\n\n");
}


template<typename MatrixType>
void print_crs_matrix(const MatrixType& matrix, const char* label = "CRS Matrix") {
    
    int numRows = matrix.numRows();
    int numCols = matrix.numCols();
    
    printf("\n=== %s (%d x %d) ===\n", label, numRows, numCols);
    
    // Print column headers
    printf("    ");
    for (int j = 0; j < numCols; j++) {
        printf("%9d ", j);
    }
    printf("\n");
    
    // Simple host loop - let Kokkos handle device access
    for (int i = 0; i < numRows; i++) {
        printf("%2d: ", i);

        // Get sparse row view for this row
        auto row = matrix.row(i);
        
        // For each column, search for the value in sparse row
        for (int j = 0; j < numCols; j++) {
            float value = 0.0f;  // Default to zero
            
            // Search through this row's entries
            for (int k = 0; k < row.length; k++) {
                if (row.colidx(k) == j) {
                    value = row.value(k);
                    break;
                }
            }
            
            printf("%9.6f ", value);
        }
        printf("\n");
    }
    
    printf("========================\n\n");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxFFKokkos<DeviceType>::calculate_H_k(const double &r, const double &shld) const
{
  double taper = tap7;
  taper = Kokkos::fma(r, taper, tap6);
  taper = Kokkos::fma(r, taper, tap5);
  taper = Kokkos::fma(r, taper, tap4);
  taper = Kokkos::fma(r, taper, tap3);
  taper = Kokkos::fma(r, taper, tap2);
  taper = Kokkos::fma(r, taper, tap1);
  taper = Kokkos::fma(r, taper, tap0);
  return taper * EV_TO_KCAL_PER_MOL / cbrt(r*r*r + shld);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::build_crs_matrix()
{
  if (crs_matrix_allocated) {
    crs_matrix = CRSMatrixType();
  }

  int nlocal = atomKK->nlocal;
  auto d_tag = atomKK->k_tag.template view<DeviceType>();
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  constexpr double EPSILON = 0.0001;

  // Count phase - store counts per row
  Kokkos::View<int*, DeviceType> d_row_nnz("row_nnz", nlocal);
  
  Kokkos::parallel_for("Count", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        int count = 1; // diagonal
        const int jnum = d_numneigh[i];
        
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj) & NEIGHMASK;
          const double dx = d_x(j,0) - d_x(i,0);
          const double dy = d_x(j,1) - d_x(i,1);
          const double dz = d_x(j,2) - d_x(i,2);
          const double rsq = dx*dx + dy*dy + dz*dz;
          
          if (rsq <= cutsq) {
            bool include = (j < nlocal) || 
                         (d_tag[i] < d_tag[j]) ||
                         (d_tag[i] == d_tag[j] && 
                          (dz > EPSILON || 
                           (fabs(dz) < EPSILON && 
                            (dy > EPSILON || 
                             (fabs(dy) < EPSILON && dx > EPSILON)))));
            if (include) count++;
          }
        }
        d_row_nnz(i) = count;
      }
    });
  
  // Build proper row map with exclusive scan
  typename CRSMatrixType::row_map_type::non_const_type row_map("crs_row_map", nlocal + 1);
  
  Kokkos::parallel_for("InitRowMap", 1, KOKKOS_LAMBDA(const int) {
    row_map(0) = 0;
  });
  
  Kokkos::parallel_scan("BuildRowMap", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i, typename CRSMatrixType::size_type& partial_sum, const bool final) {
      if (final) row_map(i+1) = partial_sum + d_row_nnz(i);
      partial_sum += d_row_nnz(i);
    });
  
  typename CRSMatrixType::size_type total_nnz;
  Kokkos::deep_copy(total_nnz, Kokkos::subview(row_map, nlocal));
  
  // Allocate final storage
  typename CRSMatrixType::index_type::non_const_type columns("crs_columns", total_nnz);
  typename CRSMatrixType::values_type::non_const_type values("crs_values", total_nnz);

  // Fill phase - each thread fills its assigned rows
  Kokkos::parallel_for("Fill", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        const int itype = d_type(i);
        int pos = row_map(i);
        
        // Diagonal first
        columns(pos) = i;
        values(pos) = d_params(itype).eta;
        pos++;
        
        // Off-diagonal
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj) & NEIGHMASK;
          const double dx = d_x(j,0) - d_x(i,0);
          const double dy = d_x(j,1) - d_x(i,1);
          const double dz = d_x(j,2) - d_x(i,2);
          const double rsq = dx*dx + dy*dy + dz*dz;
          
          if (rsq <= cutsq) {
            bool include = (j < nlocal) || 
                         (d_tag[i] < d_tag[j]) ||
                         (d_tag[i] == d_tag[j] && 
                          (dz > EPSILON || 
                           (fabs(dz) < EPSILON && 
                            (dy > EPSILON || 
                             (fabs(dy) < EPSILON && dx > EPSILON)))));
            
            if (include) {
              const double r = sqrt(rsq);
              const double shldij = d_shield(itype, d_type(j));
              columns(pos) = j;
              values(pos) = calculate_H_k(r, shldij);
              pos++;
            }
          }
        }
      }
    });
  
  crs_matrix = CRSMatrixType("H_crs", nlocal, nlocal, total_nnz, values, row_map, columns);
  crs_matrix_allocated = true;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::update_crs_matrix_values()
{
  // Update only the matrix values without changing sparsity pattern
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();
  
  int nlocal = atomKK->nlocal;
  auto values = crs_matrix.values;
  auto row_map = crs_matrix.graph.row_map;
  auto columns = crs_matrix.graph.entries;

  // Update all matrix entries
  Kokkos::parallel_for("UpdateCRSValues", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      const int itype = d_type(i);
      
      // Get row bounds
      const int row_start = row_map(i);
      const int row_end = row_map(i+1);
      
      // Process each entry in the row
      for (int pos = row_start; pos < row_end; pos++) {
        const int j = columns(pos);
        
        if (j == i) {
          // Diagonal element
          values(pos) = d_params(itype).eta;
        } else {
          // Off-diagonal element - recalculate H value
          const double delx = d_x(j,0) - d_x(i,0);
          const double dely = d_x(j,1) - d_x(i,1);
          const double delz = d_x(j,2) - d_x(i,2);
          const double rsq = delx*delx + dely*dely + delz*delz;
          
          const double r = sqrt(rsq);
          const double shldij = d_shield(itype, d_type(j));
          values(pos) = calculate_H_k(r, shldij);
        }
      }
    });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init_matvec()
{
  auto d_q = atomKK->k_q.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.view<DeviceType>();

  Kokkos::parallel_for("init_matvec", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      const int itype = d_type(i);
      if (d_mask[i] & groupbit) {
        d_Hdia_inv(i) = 1.0 / d_params(itype).eta;

        // Two separate RHS vectors
        d_rhs(i, 0) = -d_params(itype).chi - d_chi_field[i];
        d_rhs(i, 1) = -1.0;
        
        // Always use theta for extended Lagrangian initialization
        if (update->ntimestep > 0) {
          // Use theta to approximate initial s and t
          // Since q = s - u*t and theta ≈ q, use simple initialization
          d_sol(i, 0) = d_theta(i);
          d_sol(i, 1) = 1.0;  // Neutral guess
        } else {
          // Standard initialization
          d_sol(i, 0) = d_q(i);
          d_sol(i, 1) = 0.0;
        }
      }
    });
}

/* ---------------------------------------------------------------------- */

// Move this OUTSIDE of cg_solve() function - place it in the class scope
// Use the correct template signature that KokkosBatched CG expects
template<class DeviceType>
struct BatchedMatVec {

    typedef float compute_t;
    // Use the actual LAMMPS types - adjust these to match your class members
    using CRSMatrixType = KokkosSparse::CrsMatrix<compute_t, int, DeviceType>;  // or your matrix type
    using t_compute_2d = Kokkos::View<compute_t**, typename DeviceType::array_layout, DeviceType>;
    
    CRSMatrixType matrix;
    int nlocal, nall;

    BatchedMatVec(const CRSMatrixType& A, int nl, int na) 
      : matrix(A), nlocal(nl), nall(na) {}
    
    // Method 1: Basic apply (3 arguments) - output = A * input
    // Use int template parameters instead of class types to avoid C++20 requirement
    template<int TransType, int ModeType, typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void apply(const MemberType& member, 
               const t_compute_2d& input, 
               t_compute_2d& output) const {
        
        // Use team-level parallelism for matrix-vector multiply
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, nlocal), [&](const int i) {
            const int row_start = matrix.graph.row_map(i);
            const int row_end = matrix.graph.row_map(i+1);
            
            compute_t sum0 = 0.0;
            compute_t sum1 = 0.0;
            
            for (int k = row_start; k < row_end; k++) {
                const int j = matrix.graph.entries(k);
                if (j < nall) {
                    sum0 = Kokkos::fma(matrix.values(k), input(j, 0), sum0);
                    sum1 = Kokkos::fma(matrix.values(k), input(j, 1), sum1);
                }
            }
            
            output(i, 0) = sum0;
            output(i, 1) = sum1;
        });
        
        // Synchronize team members
        member.team_barrier();
    }
    
    // Method 2: AXPY apply (5 arguments) - output = alpha * A * input + beta * output
    template<int TransType, int ModeType, typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void apply(const MemberType& member, 
               const t_compute_2d& input, 
               t_compute_2d& output,
               const compute_t alpha,
               const compute_t beta) const {
        
        // Use team-level parallelism for matrix-vector multiply with AXPY
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, nlocal), [&](const int i) {
            const int row_start = matrix.graph.row_map(i);
            const int row_end = matrix.graph.row_map(i+1);
            
            compute_t sum0 = 0.0;
            compute_t sum1 = 0.0;
            
            for (int k = row_start; k < row_end; k++) {
                const int j = matrix.graph.entries(k);
                if (j < nall) {
                    sum0 = Kokkos::fma(matrix.values(k), input(j, 0), sum0);
                    sum1 = Kokkos::fma(matrix.values(k), input(j, 1), sum1);
                }
            }
            
            // AXPY operation: output = alpha * A * input + beta * output
            output(i, 0) = alpha * sum0 + beta * output(i, 0);
            output(i, 1) = alpha * sum1 + beta * output(i, 1);
        });
        
        // Synchronize team members
        member.team_barrier();
    }
};

// ALTERNATIVE APPROACH: If the template approach still fails, try this simpler version
// that accepts any template arguments but ignores them:
/*
template<class DeviceType>
struct SimpleBatchedMatVec {
    CRSMatrixType matrix;
    int nlocal, nall;

    SimpleBatchedMatVec(const CRSMatrixType& A, int nl, int na) 
      : matrix(A), nlocal(nl), nall(na) {}
    
    // Catch-all template method that accepts any template arguments
    template<typename... Args, typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void apply(const MemberType& member, 
               const t_compute_2d& input, 
               t_compute_2d& output) const {
        // Your 3-argument matrix-vector multiply here
    }
    
    template<typename... Args, typename MemberType>
    KOKKOS_INLINE_FUNCTION  
    void apply(const MemberType& member,
               const t_compute_2d& input,
               t_compute_2d& output,
               const compute_t alpha,
               const compute_t beta) const {
        // Your 5-argument AXPY matrix-vector multiply here  
    }
};
*/

// Then in your cg_solve() function, use it like this:
// BatchedMatVec<DeviceType> matvec_op(crs_matrix, nlocal, nall);



// Then in your cg_solve() function, use it like this:
// BatchedMatVec<DeviceType> matvec_op(crs_matrix, nlocal, nall);

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::cg_solve()
{
  int nlocal = atomKK->nlocal;
  int nall = nlocal + atomKK->nghost;

  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  // Custom preconditioner functor
  struct DiagonalPreconditioner {
    t_compute_1d diag_inv;
    int nn_local;
    
    DiagonalPreconditioner(t_compute_1d &d, int n) : diag_inv(d), nn_local(n) {}

    KOKKOS_FUNCTION
    void operator()(const t_compute_2d& r, t_compute_2d& z) const {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn_local), [&](const int i) {
        for (int col = 0; col < r.extent(1); col++) {
          z(i, col) = r(i, col) * diag_inv(i);
        }
      });
    }
  };
  
  // Set up KokkosBatched::CG solver
  BatchedMatVec<DeviceType> matvec_op(crs_matrix, nlocal, nall);
  DiagonalPreconditioner precond_op(d_Hdia_inv, nn);

  // CG solver parameters

  using Layout = typename DeviceType::array_layout;
  using NormViewType = Kokkos::View<accumulate_t**, Layout, DeviceType>;
  using IntViewType = Kokkos::View<int*, Layout, DeviceType>;
  using ViewType3D = Kokkos::View<compute_t***, Layout, DeviceType>;
  using KrylovHandleType = KokkosBatched::KrylovHandle<NormViewType, IntViewType, ViewType3D>;
  KrylovHandleType handle(2, 1, maxiter, false);
  handle.set_max_iteration(maxiter);
  handle.set_tolerance(tolerance);

  using ExecutionSpace = typename DeviceType::execution_space;
  using policy_type = Kokkos::TeamPolicy<ExecutionSpace>;
  policy_type policy(2, Kokkos::AUTO());

  // Use KOKKOS_LAMBDA instead of functor - much cleaner!
  Kokkos::parallel_for("CG_Solve", policy, KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

    KokkosBatched::TeamCG<typename policy_type::member_type>::invoke(
        member,
        matvec_op,
        d_rhs,
        d_sol,
        handle
    );
  });

  // FIXME: preconditioner
  // FIXME: int total_iter = result.num_iterations;


/*
 // Handle communication if needed
  if (neighflag != FULL) {
    // Convert back to separate views for communication
    Kokkos::parallel_for("extract_s", Kokkos::RangePolicy<DeviceType>(0, nall),
      KOKKOS_LAMBDA(const int i) {
        d_s(i) = static_cast<double>(d_sol(i, 0));
      });
    pack_flag = 1;
    comm->forward_comm(this);
    
    Kokkos::parallel_for("extract_t", Kokkos::RangePolicy<DeviceType>(0, nall),
      KOKKOS_LAMBDA(const int i) {
        d_t(i) = static_cast<double>(d_sol(i, 1));
      });
    pack_flag = 2;
    comm->forward_comm(this);
  }

*/


  
  // Convert final solution back to double precision in d_s and d_t
  Kokkos::parallel_for("finalize_batched_cg", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        d_s(i) = static_cast<double>(d_sol(i, 0));
        d_t(i) = static_cast<double>(d_sol(i, 1));
      }
    });

  // FIXME: compute_scalar()
  //return total_iter;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
  auto d_q = atomKK->k_q.view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  // Calculate local sums using Kokkos reductions
  double s_sum_local = 0.0, t_sum_local = 0.0;
  
  Kokkos::parallel_reduce("sum_s_t", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii, double& s_local, double& t_local) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        s_local += d_s(i);
        t_local += d_t(i);
      }
    }, s_sum_local, t_sum_local);
  
  // Calculate local u parameter with target charge consideration
  double u_local = (fabs(target_charge) > QSUMSMALL) ? 
                   (s_sum_local - target_charge) / t_sum_local :
                   s_sum_local / t_sum_local;

  // Calculate charges locally: q = s - u*t
  Kokkos::parallel_for("calculate_charges", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        d_q(i) = Kokkos::fma(-u_local, d_t(i), d_s(i));
      }
    });

  // Single MPI reduction to synchronize charges across all processes
  double* q_ptr = d_q.data();
  MPI_Allreduce(MPI_IN_PLACE, q_ptr, atom->nlocal, MPI_DOUBLE, MPI_SUM, world);

  atomKK->modified(execution_space, Q_MASK);
  
  // Forward communicate charges
  pack_flag = 4;
  comm->forward_comm(this);
}
