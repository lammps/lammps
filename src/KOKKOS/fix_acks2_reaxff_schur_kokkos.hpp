// clang-format off
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

#include "KokkosSparse_BsrMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas1_dot.hpp"
#include "KokkosBlas1_axpby.hpp"


template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::allocate_bsr_matrix()
{
  // Count total number of blocks needed
  typename AT::t_int_scalar total_blocks("total_blocks");
  typename AT::t_int_scalar h_blocks("h_blocks");
  typename AT::t_int_scalar x_blocks("x_blocks");

  Kokkos::View<int*, DeviceType> block_counts("block_counts", nn);
  
  // Count blocks from neighbor interactions
  Kokkos::parallel_for("acks2/bsr:count_blocks", nn, KOKKOS_LAMBDA(const int ii) {
    int h_count = 0;
    int x_count = 0;
    const int i = d_ilist[ii];
    
    if (d_mask[i] & groupbit) {
      const int itype = d_type(i);
      const F_FLOAT xtmp = d_x(i,0);
      const F_FLOAT ytmp = d_x(i,1);
      const F_FLOAT ztmp = d_x(i,2);
      const tagint itag = d_tag(i);
      const int jnum = d_numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = d_neighbors(i,jj);
        j &= NEIGHMASK;
        
        const F_FLOAT delx = d_x(j,0) - xtmp;
        const F_FLOAT dely = d_x(j,1) - ytmp;
        const F_FLOAT delz = d_x(j,2) - ztmp;
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        
        if (rsq <= cutsq) {
          h_count++;
          
          const int jtype = d_type(j);
          const F_FLOAT bcutoff = d_bcut(itype,jtype);
          if (rsq <= bcutoff * bcutoff) {
            x_count++;
          }
        }
      }
    }
    
    // Each neighbor interaction creates one 2x2 block
    block_counts(ii) = (h_count > x_count) ? h_count : x_count;
  });

  // Sum up total blocks
  int nnz_blocks = 0;
  Kokkos::parallel_reduce("acks2/bsr:total_blocks", nn, KOKKOS_LAMBDA(const int ii, int& sum) {
    sum += block_counts(ii);
  }, nnz_blocks);
  
  // Add diagonal blocks
  nnz_blocks += nlocal;
  
  // Add atom-global interaction blocks  
  nnz_blocks += nlocal;
  
  // Add global row blocks
  if (last_rows_flag) {
    nnz_blocks += nlocal + 1; // Interactions with all atoms + diagonal
  }
  
  // Allocate BSR matrix structure with safety margin
  nnz_blocks = static_cast<int>(nnz_blocks * 1.1);
  
  int num_block_rows = nlocal + 1; // +1 for global row
  d_row_block_counts = typename AT::t_int_1d("row_block_counts", num_block_rows);
  
  // Don't create the BSR matrix here - we'll create it after filling the data
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::build_bsr_matrix()
{
  int num_block_rows = nlocal + 1; // +1 for global row
  
  // First pass: count blocks per row
  Kokkos::deep_copy(d_row_block_counts, 0);
  
  // Count diagonal and atom-global blocks
  Kokkos::parallel_for("acks2/bsr:count_base_blocks", nlocal, KOKKOS_LAMBDA(const int i) {
    Kokkos::atomic_add(&d_row_block_counts(i), 2); // diagonal + atom-global
  });
  
  // Count neighbor blocks
  Kokkos::parallel_for("acks2/bsr:count_neighbor_blocks", nn, KOKKOS_LAMBDA(const int ii) {
    const int i = d_ilist[ii];
    
    if (d_mask[i] & groupbit) {
      const int itype = d_type(i);
      const F_FLOAT xtmp = d_x(i,0);
      const F_FLOAT ytmp = d_x(i,1);
      const F_FLOAT ztmp = d_x(i,2);
      const tagint itag = d_tag(i);
      const int jnum = d_numneigh[i];
      
      int neighbor_blocks = 0;
      
      for (int jj = 0; jj < jnum; jj++) {
        int j = d_neighbors(i,jj);
        j &= NEIGHMASK;
        
        if (neighflag != FULL) {
          const tagint jtag = d_tag(j);
          if (j >= nlocal) {
            if (itag > jtag) {
              if ((itag+jtag) % 2 == 0) continue;
            } else if (itag < jtag) {
              if ((itag+jtag) % 2 == 1) continue;
            } else {
              if (d_x(j,2) < ztmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) < ytmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) == ytmp && d_x(j,0) < xtmp) continue;
            }
          }
        }
        
        const F_FLOAT delx = d_x(j,0) - xtmp;
        const F_FLOAT dely = d_x(j,1) - ytmp;
        const F_FLOAT delz = d_x(j,2) - ztmp;
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        
        if (rsq <= cutsq) {
          neighbor_blocks++;
        }
      }
      
      Kokkos::atomic_add(&d_row_block_counts(i), neighbor_blocks);
    }
  });
  
  // Count global row blocks
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/bsr:count_global_blocks", 1, KOKKOS_LAMBDA(const int) {
      d_row_block_counts(nlocal) = nlocal + 1; // All atoms + diagonal
    });
  }
  
  // Second pass: compute block_row_map
  // Create non-const views that we'll fill
  typename BSRMatrixType::row_map_type::non_const_type row_map("block_row_map", num_block_rows + 1);
  
  Kokkos::parallel_scan("acks2/bsr:block_row_map", num_block_rows,
    KOKKOS_LAMBDA(const int i, int& update, const bool final) {
      const int count = d_row_block_counts(i);
      if (final) {
        row_map(i) = update;
      }
      update += count;
    });
  
  // Get total blocks
  typename AT::t_int_scalar d_total_blocks("total_blocks");
  Kokkos::deep_copy(d_total_blocks, Kokkos::subview(row_map, num_block_rows));
  int h_total_blocks;
  Kokkos::deep_copy(h_total_blocks, d_total_blocks);
  
  Kokkos::parallel_for("acks2/bsr:set_final_row_map", 1, KOKKOS_LAMBDA(const int) {
    row_map(num_block_rows) = h_total_blocks;
  });
  
  // Now allocate cols and values with the correct size
  typename BSRMatrixType::index_type::non_const_type cols("block_cols", h_total_blocks);
  typename BSRMatrixType::values_type::non_const_type values("values", h_total_blocks * 4);
  
  // Third pass: fill matrix data
  Kokkos::View<int*, DeviceType> d_block_pos("block_positions", num_block_rows);
  Kokkos::deep_copy(d_block_pos, row_map);
  
  // Fill diagonal blocks  
  Kokkos::parallel_for("acks2/bsr:diagonal_blocks", nlocal, KOKKOS_LAMBDA(const int i) {
    const int itype = d_type(i);
    const double eta = params(itype).eta;
    const double x_diag = d_X_diag[i];
    
    int pos = Kokkos::atomic_fetch_add(&d_block_pos(i), 1);
    int value_offset = pos * 4; // Each 2x2 block has 4 values
    
    cols(pos) = i;
    values(value_offset + 0) = eta;      // (0,0)
    values(value_offset + 1) = 1.0;      // (0,1)
    values(value_offset + 2) = 1.0;      // (1,0)  
    values(value_offset + 3) = x_diag;   // (1,1)
  });
  
  // Fill neighbor blocks
  Kokkos::parallel_for("acks2/bsr:neighbor_blocks", nn, KOKKOS_LAMBDA(const int ii) {
    const int i = d_ilist[ii];
    
    if (d_mask[i] & groupbit) {
      const int itype = d_type(i);
      const F_FLOAT xtmp = d_x(i,0);
      const F_FLOAT ytmp = d_x(i,1);
      const F_FLOAT ztmp = d_x(i,2);
      const tagint itag = d_tag(i);
      const int jnum = d_numneigh[i];
      
      for (int jj = 0; jj < jnum; jj++) {
        int j = d_neighbors(i,jj);
        j &= NEIGHMASK;
        const int jtype = d_type(j);
        
        if (neighflag != FULL) {
          const tagint jtag = d_tag(j);
          if (j >= nlocal) {
            if (itag > jtag) {
              if ((itag+jtag) % 2 == 0) continue;
            } else if (itag < jtag) {
              if ((itag+jtag) % 2 == 1) continue;
            } else {
              if (d_x(j,2) < ztmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) < ytmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) == ytmp && d_x(j,0) < xtmp) continue;
            }
          }
        }
        
        const F_FLOAT delx = d_x(j,0) - xtmp;
        const F_FLOAT dely = d_x(j,1) - ytmp;
        const F_FLOAT delz = d_x(j,2) - ztmp;
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        
        if (rsq <= cutsq) {
          const F_FLOAT r = sqrt(rsq);
          const F_FLOAT shldij = d_shield(itype,jtype);
          double h_val = calculate_H_k(r, shldij);
          
          double x_val = 0.0;
          const F_FLOAT bcutoff = d_bcut(itype,jtype);
          if (rsq <= bcutoff * bcutoff) {
            x_val = calculate_X_k(r, bcutoff);
          }
          
          int pos = Kokkos::atomic_fetch_add(&d_block_pos(i), 1);
          int value_offset = pos * 4;
          
          cols(pos) = j;
          values(value_offset + 0) = h_val;    // (0,0)
          values(value_offset + 1) = 0.0;      // (0,1)
          values(value_offset + 2) = 0.0;      // (1,0)
          values(value_offset + 3) = x_val;    // (1,1)
        }
      }
    }
  });
  
  // Fill atom-global interaction blocks
  Kokkos::parallel_for("acks2/bsr:atom_global_blocks", nlocal, KOKKOS_LAMBDA(const int i) {
    int pos = Kokkos::atomic_fetch_add(&d_block_pos(i), 1);
    int value_offset = pos * 4;
    
    cols(pos) = nlocal; // Global block index
    values(value_offset + 0) = 0.0;      // (0,0)
    values(value_offset + 1) = 1.0;      // (0,1) s_H - λ_2
    values(value_offset + 2) = 1.0;      // (1,0) s_X - λ_1
    values(value_offset + 3) = 0.0;      // (1,1)
  });
  
  // Fill global row blocks
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/bsr:global_row_blocks", nlocal, KOKKOS_LAMBDA(const int i) {
      int pos = Kokkos::atomic_fetch_add(&d_block_pos(nlocal), 1);
      int value_offset = pos * 4;
      
      cols(pos) = i;
      values(value_offset + 0) = 0.0;      // (0,0)
      values(value_offset + 1) = 1.0;      // (0,1)
      values(value_offset + 2) = 1.0;      // (1,0)
      values(value_offset + 3) = 0.0;      // (1,1)
    });
    
    // Global diagonal block
    Kokkos::parallel_for("acks2/bsr:global_diagonal", 1, KOKKOS_LAMBDA(const int) {
      int pos = Kokkos::atomic_fetch_add(&d_block_pos(nlocal), 1);
      int value_offset = pos * 4;
      
      cols(pos) = nlocal;
      values(value_offset + 0) = 0.0;      // (0,0)
      values(value_offset + 1) = 0.0;      // (0,1)
      values(value_offset + 2) = 0.0;      // (1,0)
      values(value_offset + 3) = 0.0;      // (1,1)
    });
  }
  
  // Create the BSR matrix with the filled data
  bsr_matrix = BSRMatrixType("A_bsr", num_block_rows, num_block_rows, h_total_blocks,
                            values, row_map, cols, 2);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::schur_cg_solve()
{
  const int n = nlocal;
  const int block_size = 2;
  const int maxiter = imax;
  const double tol = tolerance;
  
  // For mixed precision optimization, we'll use float for some intermediate calculations
  typename AT::t_float_1d d_r2("r2", 2);      // Float precision residual
  typename AT::t_float_1d d_p2("p2", 2);      // Float precision search direction
  typename AT::t_float_1d d_Ap2("Ap2", 2);    // Float precision A*p
  
  // Double precision for final solution and key operations
  typename AT::t_float_1d d_lambda("lambda", 2);
  typename AT::t_float_1d d_r("r", 2);
  typename AT::t_float_1d d_p("p", 2);
  typename AT::t_float_1d d_Ap("Ap", 2);
  typename AT::t_float_1d d_temp1("temp1", n * block_size);
  typename AT::t_float_1d d_temp2("temp2", n * block_size);
  
  // Extract subviews for block structure
  auto d_b1 = Kokkos::subview(d_b_s, Kokkos::make_pair(0, n * block_size));
  auto d_b2 = Kokkos::subview(d_b_s, Kokkos::make_pair(n * block_size, n * block_size + 2));
  
  // Step 1: Compute A^{-1}b₁ using diagonal preconditioner
  apply_block_diagonal_preconditioner(d_b1, d_temp1);
  
  // Step 2: Form Schur complement RHS: r = b₂ - C*A^{-1}*b₁
  compute_C_times_vector(d_temp1, d_r);
  
  // Need to factor in target charge here
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/schur:adjust_rhs", 1, KOKKOS_LAMBDA(const int) {
      d_r(0) = target_charge - d_r(0);    // λ₁ constraint (total X)
      d_r(1) = target_charge - d_r(1);    // λ₂ constraint (total H = total charge)
    });
  } else {
    Kokkos::deep_copy(d_r, 0.0);
  }
  
  // Initialize CG iteration
  Kokkos::deep_copy(d_lambda, 0.0);
  Kokkos::deep_copy(d_p, d_r);
  
  double rho = KokkosBlas::dot(d_r, d_r);
  double rho0 = rho;
  
  int iter = 0;
  
  while (rho > tol * tol * rho0 && iter < maxiter) {
    // Compute Ap = S*p = -C*A^{-1}*B*p
    compute_B_times_vector(d_p, d_temp2);
    apply_block_diagonal_preconditioner(d_temp2, d_temp1);
    compute_C_times_vector(d_temp1, d_Ap);
    
    // Schur complement is negative of this
    Kokkos::parallel_for("acks2/schur:negate_Ap", 2, KOKKOS_LAMBDA(const int i) {
      d_Ap(i) = -d_Ap(i);
    });
    
    double pAp = KokkosBlas::dot(d_p, d_Ap);
    double alpha = rho / pAp;
    
    // Update solution: λ = λ + alpha*p
    KokkosBlas::axpy(alpha, d_p, d_lambda);
    
    // Update residual: r = r - alpha*Ap
    KokkosBlas::axpy(-alpha, d_Ap, d_r);
    
    double rho_new = KokkosBlas::dot(d_r, d_r);
    double beta = rho_new / rho;
    
    // Update search direction: p = r + beta*p
    Kokkos::parallel_for("acks2/schur:update_p", 2, KOKKOS_LAMBDA(const int i) {
      d_p(i) = d_r(i) + beta * d_p(i);
    });
    
    rho = rho_new;
    iter++;
  }
  
  // Step 4: Recover solution x₁ = A^{-1}(b₁ - B*λ)
  compute_B_times_vector(d_lambda, d_temp2);
  
  Kokkos::parallel_for("acks2/schur:recover_solution", n * block_size, KOKKOS_LAMBDA(const int i) {
    d_temp2(i) = d_b1(i) - d_temp2(i);
  });
  
  apply_block_diagonal_preconditioner(d_temp2, d_temp1);
  
  // Copy solution back
  Kokkos::parallel_for("acks2/schur:copy_solution", n, KOKKOS_LAMBDA(const int i) {
    d_s(i) = d_temp1(i * 2);            // s_H component
    d_s(n + i) = d_temp1(i * 2 + 1);    // s_X component
  });
  
  // Copy global variables
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/schur:copy_global", 1, KOKKOS_LAMBDA(const int) {
      d_s(2*n) = d_lambda(0);      // λ₁
      d_s(2*n + 1) = d_lambda(1);  // λ₂
    });
  }
  
  if (comm->me == 0 && iter >= maxiter) {
    error->warning(FLERR, "Fix acks2/reaxff/kk Schur CG convergence failed after {} iterations", iter);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::apply_block_diagonal_preconditioner(
    const typename AT::t_float_1d& d_in,
    typename AT::t_float_1d& d_out)
{
  const int n = nlocal;
  
  Kokkos::parallel_for("acks2/precond:block_diagonal", n, KOKKOS_LAMBDA(const int i) {
    const int itype = d_type(i);
    const double h_diag = params(itype).eta;
    const double x_diag = d_X_diag[i];
    
    // 2x2 block inverse: [eta, 1; 1, x_diag]^{-1}
    const double det = h_diag * x_diag - 1.0;
    
    const double in1 = d_in(i * 2);
    const double in2 = d_in(i * 2 + 1);
    
    // Apply inverse
    d_out(i * 2) = (x_diag * in1 - in2) / det;
    d_out(i * 2 + 1) = (-in1 + h_diag * in2) / det;
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::compute_B_times_vector(
    const typename AT::t_float_1d& d_in,
    typename AT::t_float_1d& d_out)
{
  const int n = nlocal;
  
  // B = [ 0  I ]
  //     [ I  0 ]
  
  Kokkos::parallel_for("acks2/B_times_vec", n, KOKKOS_LAMBDA(const int i) {
    d_out(i * 2) = d_in(1);        // s_H couples to λ₂
    d_out(i * 2 + 1) = d_in(0);    // s_X couples to λ₁
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::compute_C_times_vector(
    const typename AT::t_float_1d& d_in,
    typename AT::t_float_1d& d_out)
{
  const int n = nlocal;
  
  // C = [ 0  I ]
  //     [ I  0 ]
  
  typename AT::t_float_scalar sum0("sum0");
  typename AT::t_float_scalar sum1("sum1");
  
  Kokkos::parallel_reduce("acks2/C_times_vec", n, KOKKOS_LAMBDA(const int i, double& lsum0, double& lsum1) {
    lsum0 += d_in(i * 2 + 1);  // Sum of s_X
    lsum1 += d_in(i * 2);      // Sum of s_H
  }, sum0, sum1);
  
  // MPI reduction
  double global_sum0, global_sum1;
  Kokkos::deep_copy(global_sum0, sum0);
  Kokkos::deep_copy(global_sum1, sum1);
  
  MPI_Allreduce(MPI_IN_PLACE, &global_sum0, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &global_sum1, 1, MPI_DOUBLE, MPI_SUM, world);
  
  // Set output
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/C_set_output", 1, KOKKOS_LAMBDA(const int) {
      d_out(0) = global_sum0;
      d_out(1) = global_sum1;
    });
  } else {
    Kokkos::deep_copy(d_out, 0.0);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::sparse_matvec_bsr(
    const typename AT::t_float_1d& x_in,
    typename AT::t_float_1d& b_out)
{
  KokkosSparse::spmv("N", 1.0, bsr_matrix, x_in, 0.0, b_out);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  atomKK->sync(execution_space, datamask_read);

  d_x = atomKK->k_x.view<DeviceType>();
  d_q = atomKK->k_q.view<DeviceType>();
  d_tag = atomKK->k_tag.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  nlocal = atomKK->nlocal;
  newton_pair = force->newton_pair;

  k_params.template sync<DeviceType>();
  k_shield.template sync<DeviceType>();
  k_bcut.template sync<DeviceType>();
  k_tap.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  nn = list->inum;
  NN = atom->nlocal + atom->nghost;

  copymode = 1;

  // Allocate arrays if needed
  allocate_array();

  // Check if we need to rebuild matrix
  if (!allocated_flag || last_allocate < neighbor->lastcall) {
    allocate_bsr_matrix();
    last_allocate = update->ntimestep;
    
    // Set up last_rows_flag
    int flag = comm->me;
    if (nn == 0) flag = MAXSMALLINT;
    MPI_Allreduce(&flag, &last_rows_rank, 1, MPI_INT, MPI_MIN, world);
    last_rows_flag = (comm->me == last_rows_rank);
  }
  
  // Compute X diagonal elements  
  compute_X_diagonal();
  
  // Build BSR matrix
  build_bsr_matrix();
  
  // Update chi field if needed
  if (efield) get_chi_field();
  
  // Initialize right-hand side
  init_matvec();
  
  // Solve using Schur complement CG
  schur_cg_solve();
  
  // Calculate atomic charges
  calculate_Q();
  
  copymode = 0;
  
  if (!allocated_flag) allocated_flag = 1;
  
  atomKK->modified(execution_space, datamask_modify);
  k_s.modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::compute_X_diagonal()
{
  Kokkos::deep_copy(d_X_diag, 0.0);
  
  Kokkos::parallel_for("acks2/X_diagonal", nn, KOKKOS_LAMBDA(const int ii) {
    const int i = d_ilist[ii];
    
    if (d_mask[i] & groupbit) {
      const int itype = d_type(i);
      const F_FLOAT xtmp = d_x(i,0);
      const F_FLOAT ytmp = d_x(i,1);
      const F_FLOAT ztmp = d_x(i,2);
      const tagint itag = d_tag(i);
      const int jnum = d_numneigh[i];
      
      double x_diag_sum = 0.0;
      
      for (int jj = 0; jj < jnum; jj++) {
        int j = d_neighbors(i,jj);
        j &= NEIGHMASK;
        const int jtype = d_type(j);
        
        const F_FLOAT delx = d_x(j,0) - xtmp;
        const F_FLOAT dely = d_x(j,1) - ytmp;
        const F_FLOAT delz = d_x(j,2) - ztmp;
        const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        
        if (rsq > cutsq) continue;
        
        const F_FLOAT bcutoff = d_bcut(itype,jtype);
        if (rsq > bcutoff * bcutoff) continue;
        
        if (neighflag != FULL) {
          const tagint jtag = d_tag(j);
          if (j >= nlocal) {
            if (itag > jtag) {
              if ((itag+jtag) % 2 == 0) continue;
            } else if (itag < jtag) {
              if ((itag+jtag) % 2 == 1) continue;
            } else {
              if (d_x(j,2) < ztmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) < ytmp) continue;
              if (d_x(j,2) == ztmp && d_x(j,1) == ytmp && d_x(j,0) < xtmp) continue;
            }
          }
        }
        
        const F_FLOAT r = sqrt(rsq);
        const double x_val = calculate_X_k(r, bcutoff);
        
        x_diag_sum -= x_val;
        
        if (neighflag != FULL && j < nlocal) {
          Kokkos::atomic_add(&d_X_diag[j], -x_val);
        }
      }
      
      Kokkos::atomic_add(&d_X_diag[i], x_diag_sum);
    }
  });
  
  // Communication for reverse neighbors if needed
  if (neighflag != FULL) {
    pack_flag = 4;
    k_X_diag.template modify<DeviceType>();
    k_X_diag.template sync<LMPHostType>();
    comm->reverse_comm(this);
    k_X_diag.template modify<LMPHostType>();
    k_X_diag.template sync<DeviceType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init_matvec()
{
  const int n = nlocal;
  
  // Initialize right-hand side
  Kokkos::parallel_for("acks2/init_rhs", nn, KOKKOS_LAMBDA(const int ii) {
    const int i = d_ilist[ii];
    const int itype = d_type(i);
    
    if (d_mask[i] & groupbit) {
      d_b_s[i] = -params(itype).chi - d_chi_field[i];
      d_b_s[n + i] = 0.0;
    }
  });
  
  // Initialize global constraints
  if (last_rows_flag) {
    Kokkos::parallel_for("acks2/init_global_rhs", 1, KOKKOS_LAMBDA(const int) {
      d_b_s[2*n] = target_charge;      // Total X = target_charge
      d_b_s[2*n + 1] = target_charge;  // Total H = target_charge
    });
  }
  
  // Initialize solution vector to zero
  Kokkos::deep_copy(d_s, 0.0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::calculate_Q()
{
  const int n = nlocal;
  
  // Need to ensure solution is available on all processors
  pack_flag = 2;
  k_s.template modify<DeviceType>();
  k_s.template sync<LMPHostType>();
  comm->forward_comm(this);
  k_s.template modify<LMPHostType>();
  k_s.template sync<DeviceType>();
  
  // Extract atomic charges from solution
  Kokkos::parallel_for("acks2/extract_charges", nn, KOKKOS_LAMBDA(const int ii) {
    const int i = d_ilist[ii];
    
    if (d_mask[i] & groupbit) {
      d_q(i) = d_s(i);  // Charge is the s_H component
    }
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixACKS2ReaxFFKokkos<DeviceType>::calculate_H_k(const F_FLOAT &r, const F_FLOAT &shld) const
{
  F_FLOAT taper, denom;

  taper = d_tap[7] * r + d_tap[6];
  taper = taper * r + d_tap[5];
  taper = taper * r + d_tap[4];
  taper = taper * r + d_tap[3];
  taper = taper * r + d_tap[2];
  taper = taper * r + d_tap[1];
  taper = taper * r + d_tap[0];

  denom = r * r * r + shld;
  denom = cbrt(denom);

  return taper * EV_TO_KCAL_PER_MOL / denom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixACKS2ReaxFFKokkos<DeviceType>::calculate_X_k(const F_FLOAT &r, const F_FLOAT &bcut) const
{
  const F_FLOAT d = r/bcut;
  const F_FLOAT d3 = d*d*d;
  const F_FLOAT omd = 1.0 - d;
  const F_FLOAT omd2 = omd*omd;
  const F_FLOAT omd6 = omd2*omd2*omd2;

  return bond_softness*d3*omd6;
}

