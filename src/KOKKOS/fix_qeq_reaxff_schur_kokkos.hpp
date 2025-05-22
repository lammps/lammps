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
#include "KokkosSparse_CrsMatrix.hpp"

using Scalar = double;
using Ordinal = int;
using Offset = size_t;
using Layout = Kokkos::LayoutLeft;

template<typename ViewType>
void print_1d_view(const ViewType& view, const char* label, int max_elements = 10) {
    auto h_view = Kokkos::create_mirror_view(view);
    Kokkos::deep_copy(h_view, view);
    
    printf("%s[%d]: ", label, (int)view.extent(0));
    int n = std::min(max_elements, (int)view.extent(0));
    for (int i = 0; i < n; i++) {
        printf("%10.6f ", (double)h_view(i));
    }
    if (view.extent(0) > max_elements) printf("...");
    printf("\n");
}


template<typename MatrixType>
void print_crs_matrix(const MatrixType& matrix, const char* label = "CRS Matrix") {
    
    int numRows = matrix.numRows();
    int numCols = matrix.numCols();
    
    printf("\n=== %s (%d x %d) ===\n", label, numRows, numCols);
    
    // Print column headers
    printf("     ");
    for (int j = 0; j < numCols; j++) {
        printf("%10d ", j);
    }
    printf("\n");
    
    // Simple host loop - let Kokkos handle device access
    for (int i = 0; i < numRows; i++) {
        printf("%3d: ", i);
        
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
            
            printf("%10.6f ", value);
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
void FixQEqReaxFFKokkos<DeviceType>::sparse_matvec(t_compute_1d &in, t_compute_1d &out)
{
  int nlocal = atomKK->nlocal;
  int nall = nlocal + atomKK->nghost;
  
  // Initialize ALL atoms to zero (including ghosts)
  Kokkos::deep_copy(out, 0.0);
  
  // Apply matrix to local part
  auto in_view = Kokkos::subview(in, std::make_pair(0, nall));
  auto out_nlocal = Kokkos::subview(out, std::make_pair(0, nlocal));
  
  // Custom SpMV that can handle ghost columns
  Kokkos::parallel_for("SpMV", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      const int row_start = crs_matrix.graph.row_map(i);
      const int row_end = crs_matrix.graph.row_map(i+1);
      
      kk_compute sum = 0.0;
      for (int k = row_start; k < row_end; k++) {
        const int j = crs_matrix.graph.entries(k);
        sum += crs_matrix.values(k) * in_view(j);
      }
      out_nlocal(i) = sum;
    });
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
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::init_matvec()
{
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.view<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      const int itype = d_type(i);
      if (d_mask[i] & groupbit) {
        d_Hdia_inv[i] = 1.0 / d_params(itype).eta;
        d_b[i] = -d_params(itype).chi - d_chi_field[i];
        d_s(i) = d_theta(i);  // Use auxiliary charge as initial guess
        d_o(i) = 0.0;
        d_r(i) = 0.0;
        d_p(i) = 0.0;
        d_d(i) = 0.0;
      }
    });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::cg_solve()
{
  // r = b - H*s
  sparse_matvec(d_s, d_o);  // H*s -> o
  
  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    comm->reverse_comm(this);
    k_o.template sync<DeviceType>();
  }

  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  // Compute initial residual
  // r = b - o
  // d = r * Hdia_inv (preconditioned residual)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        d_r(i) = d_b(i) - d_o(i);
        d_d(i) = d_r(i) * d_Hdia_inv[i];
        d_p(i) = d_d(i);  // Initial conjugate direction is same as steepest descent
      }
    });
  
  // Compute initial residual norm
  double rnorm2 = KokkosBlas::dot(d_r, d_d);
  double initial_rnorm = sqrt(rnorm2);
  
  // Early termination if residual is already small
  if (initial_rnorm < tolerance) return;

  double alpha, beta, rnorm2_old;

  print_crs_matrix(crs_matrix);
  print_1d_view(d_Hdia_inv, "d_Hdia_inv", nn);
  print_1d_view(d_p, "d_p", nn);
  print_1d_view(d_r, "d_r", nn);
  print_1d_view(d_d, "d_d", nn);

  // Optimization: Allow adaptive number of iterations with convergence check
  for (int iter = 0; iter < maxiter; iter++) {
    // Compute alpha = r*d / (d*H*d)
    sparse_matvec(d_p, d_o);  // H*p -> o

    if (neighflag != FULL) {
      k_o.template modify<DeviceType>();
      comm->reverse_comm(this);
      k_o.template sync<DeviceType>();
    }

    // Compute dHd = d·o
    double dHd = KokkosBlas::dot(d_p, d_o);
    
    // Bug fix: Check for division by zero
    if (fabs(dHd) < COMPARE_TOLERANCE) break;
    
    alpha = rnorm2 / dHd;
    
    // Update solution and residual: s += alpha*p, r -= alpha*H*p
    rnorm2_old = rnorm2;
    
    // Optimization: Fuse vector operations to reduce memory access
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
      KOKKOS_LAMBDA(const int ii) {
        const int i = d_ilist[ii];
        if (d_mask[i] & groupbit) {
          d_s(i) += alpha * d_p(i);
          d_r(i) -= alpha * d_o(i);
          d_d(i) = d_r(i) * d_Hdia_inv[i];  // Preconditioned residual
        }
      });
    
    // Compute new residual norm
    rnorm2 = KokkosBlas::dot(d_r, d_d);
    
    // Check for convergence
    if (sqrt(rnorm2) / initial_rnorm < tolerance) break;

    // Compute beta = rnorm2 / rnorm2_old
    beta = rnorm2 / rnorm2_old;
    
    // Update conjugate direction: p = d + beta*p
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
      KOKKOS_LAMBDA(const int ii) {
        const int i = d_ilist[ii];
        if (d_mask[i] & groupbit) d_p(i) = d_d(i) + beta * d_p(i);
      });
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
  auto d_q = atomKK->k_q.view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  copymode = 1;

  // Update charges with the final solution
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) d_q(i) = d_s(i);  // Use s directly for charge
    });

  copymode = 0;

  atomKK->modified(execution_space, Q_MASK);

  print_1d_view(d_q, "d_q", nn);

  // Forward communicate charges
  pack_flag = 2;
  comm->forward_comm(this);
}







