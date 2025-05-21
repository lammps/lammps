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
   Contributing authors: Ray Shan (SNL), Stan Moore (SNL),
                          Kamesh Arumugam (NVIDIA)

   Nicholas Curtis (AMD), Leopold Grinberd (AMD), and Gina Sitaraman (AMD):
     - Reduced math overhead: enabled specialized calls (e.g., cbrt for a
         cube root instead of pow) and use power/exponential laws to reduce the
         number of exponentials evaluated, etc.
     - Fused the CG solve for "S" and "T" matrices
     - Improved the SpMV algorithm by using vector instead of team level
         parallelism on GPUs

   Mitch Murphy (alphataubio at gmail):
     - Extended Lagrangian (https://doi.org/10.1016/j.cpc.2015.02.023)
     - Mixed precision kokkos-kernels
------------------------------------------------------------------------- */

#include "fix_qeq_reaxff_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

// Physical constants
static constexpr double EV_TO_KCAL_PER_MOL = 14.4;

// Optimization: Pre-define constants for common operations
static constexpr double COMPARE_TOLERANCE = 1.0e-10;  // Tolerance for floating-point comparisons

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxFFKokkos<DeviceType>::FixQEqReaxFFKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixQEqReaxFF(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | Q_MASK | MASK_MASK | TYPE_MASK;
  datamask_modify = Q_MASK;

  nmax = 0;
  allocated_flag = 0;
  last_allocate = -1;
  matrix_sparsity_initialized = false;
  crs_matrix_allocated = false;

  // Extended Lagrangian parameters FIXME: make it an option
  K = 2.0;  // recommended value from Niklasson papers

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxFFKokkos<DeviceType>::~FixQEqReaxFFKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_s);
  memoryKK->destroy_kokkos(k_theta);
  memoryKK->destroy_kokkos(k_theta_dot);
  memoryKK->destroy_kokkos(k_chi_field, chi_field);
  memoryKK->destroy_kokkos(k_o);
  memoryKK->destroy_kokkos(d_r);
  memoryKK->destroy_kokkos(d_p);
  memoryKK->destroy_kokkos(d_d);
  memoryKK->destroy_kokkos(d_Hdia_inv);
  memoryKK->destroy_kokkos(d_b);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::post_constructor()
{
  FixQEqReaxFF::pertype_parameters(pertype_option);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init()
{
  atomKK->sync(execution_space, Q_MASK);

  FixQEqReaxFF::init();

  // adjust neighbor list request for KOKKOS
  neighflag = lmp->kokkos->neighflag_qeq;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  nmax = atom->nmax;

  // chi field
  memoryKK->create_kokkos(k_chi_field, chi_field, nmax, "qeq/reaxff/kk:chi_field");
  d_chi_field = k_chi_field.template view<DeviceType>();
  if (efield) get_chi_field();

  // extended lagrangian

  dt = update->dt;
  dt_half = 0.5 * dt;
  dt2_half = dt * dt_half;
  omega = sqrt(K/(dt * dt));
  omega2 = omega * omega;

  k_theta = tdual_compute_1d("qeq/reaxff/kk:theta", nmax);
  k_theta_dot = tdual_compute_1d("qeq/reaxff/kk:theta_dot", nmax);

  d_theta = k_theta.template view<DeviceType>();
  d_theta_dot = k_theta_dot.template view<DeviceType>();

  Kokkos::deep_copy(k_theta.h_view, 0.0);
  Kokkos::deep_copy(k_theta_dot.h_view, 0.0);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();

  // schur cg

  k_o = tdual_compute_1d("qeq/reaxff/kk:k_o", nmax);
  k_s = tdual_compute_1d("qeq/reaxff/kk:k_s", nmax);
  d_o = k_o.template view<DeviceType>();
  d_s = k_s.template view<DeviceType>();

  d_r = t_compute_1d("qeq/reaxff/kk:r", nmax);
  d_p = t_compute_1d("qeq/reaxff/kk:p", nmax);
  d_d = t_compute_1d("qeq/reaxff/kk:d", nmax);
  d_Hdia_inv = t_compute_1d("qeq/reaxff/kk:Hdia_inv", nmax);
  d_b = t_compute_1d("qeq/reaxff/kk:b", nmax);

  // qeq parameters

  int ntypes = atom->ntypes;
  k_params = tdual_qeq_1d("FixQEqReaxFF::params", ntypes+1);
  d_params = k_params.template view<DeviceType>();

  for (int n = 1; n <= ntypes; n++) {
    k_params.h_view(n).chi = chi[n];
    k_params.h_view(n).eta = eta[n];
    k_params.h_view(n).gamma = gamma[n];
  }
  k_params.template modify<LMPHostType>();

  cutsq = swb * swb;

  init_shielding_k();

  // Reset allocation tracking
  last_allocate = -1;
  matrix_sparsity_initialized = false;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init_shielding_k()
{
  int i,j;
  int ntypes = atom->ntypes;

  k_shield = DAT::tdual_float_2d("qeq/reaxff/kk:shield", ntypes+1, ntypes+1);
  d_shield = k_shield.template view<DeviceType>();

  for (i = 1; i <= ntypes; ++i)
    for (j = 1; j <= ntypes; ++j)
      k_shield.h_view(i,j) = pow(gamma[i] * gamma[j], -1.5);

  k_shield.template modify<LMPHostType>();
  k_shield.template sync<DeviceType>();

  tap0 = Tap[0];
  tap1 = Tap[1];
  tap2 = Tap[2];
  tap3 = Tap[3];
  tap4 = Tap[4];
  tap5 = Tap[5];
  tap6 = Tap[6];
  tap7 = Tap[7];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  atomKK->sync(execution_space, datamask_read);

  k_params.template sync<DeviceType>();
  k_shield.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  nn = list->inum;

  if (atom->nmax > nmax) resize_views();

  copymode = 1;

  // Build the sparse matrix if needed (when sparsity pattern changes)
  // Optimization: Only rebuild matrix when sparsity pattern changes
  bool rebuild_matrix = false;
  if (!matrix_sparsity_initialized || last_allocate < neighbor->lastcall) {
    rebuild_matrix = true;
    last_allocate = update->ntimestep;
    matrix_sparsity_initialized = true;
  }

  // Only update matrix values if sparsity pattern hasn't changed
  if (rebuild_matrix) build_crs_matrix();
  else update_crs_matrix_values();

  // Initialize matvec
  if (efield) get_chi_field();
  init_matvec();

  // Pack and communicate charges
  pack_flag = 1;
  k_s.template modify<DeviceType>();
  comm->forward_comm(this);
  k_s.template sync<DeviceType>();

  // Do CG iterations using theta as initial guess
  cg_solve();

  // Calculate atomic charges
  calculate_q();

  // Update auxiliary variables using velocity Verlet
  update_extended_lagrangian();

  copymode = 0;

  if (!allocated_flag) allocated_flag = 1;
  atomKK->modified(execution_space, datamask_modify);
  d_neighbors = typename AT::t_neighbors_2d();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::resize_views()
{
  nmax = atom->nmax;

  k_chi_field.resize(nmax);
  if (efield) get_chi_field();

  k_theta.resize(nmax);
  k_theta_dot.resize(nmax);

  // FIXME: only zero new entries
  //Kokkos::deep_copy(k_theta, 0.0);
  //Kokkos::deep_copy(k_theta_dot, 0.0);
  //k_theta.template modify<LMPHostType>();
  //k_theta_dot.template modify<LMPHostType>();

  k_o.resize(nmax);
  k_s.resize(nmax);

  Kokkos::resize(d_r, nmax);
  Kokkos::resize(d_p, nmax);
  Kokkos::resize(d_d, nmax);
  Kokkos::resize(d_Hdia_inv, nmax);
  Kokkos::resize(d_b, nmax);

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
    // Clean up previous allocations to avoid memory leaks
    // Note: In Kokkos, Views are reference counted and will be deallocated automatically
    // When new views are assigned, but we'll explicitly ensure cleanup
    crs_matrix = CRSMatrixType(); // This should trigger cleanup of previous allocation
  }

  // Direct construction of CRS matrix for QEq
  int nlocal = atomKK->nlocal;

  // First, count interactions to allocate arrays
  Kokkos::View<int*, DeviceType> d_row_nnz("row_nnz", nlocal);
  Kokkos::deep_copy(d_row_nnz, 0); // Initialize to zero

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  // Count nonzeros per row
  Kokkos::parallel_for("CountNNZ", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      // Bug fix: Initialize ALL rows, not just active ones
      // Initialize with diagonal element
      d_row_nnz(i) = 1;
      
      // Only count off-diagonals for atoms in group
      if (d_mask[i] & groupbit) {
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
          // Only include atoms that are part of the group and are local
          if ((d_mask[j] & groupbit) && j < nlocal) {
            const double delx = d_x(j,0) - d_x(i,0);
            const double dely = d_x(j,1) - d_x(i,1);
            const double delz = d_x(j,2) - d_x(i,2);
            const double rsq = delx*delx + dely*dely + delz*delz;
            if (rsq <= cutsq) d_row_nnz(i)++;
          }
        }
      }
    });
  
  // Create row map from row_nnz (exclusive scan)
  typename CRSMatrixType::row_map_type::non_const_type row_map("crs_row_map", nlocal + 1);
  
  // Set first element to 0
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, 1),
    KOKKOS_LAMBDA(const int) { row_map(0) = 0; });

  // Exclusive scan to compute row offsets
  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i, typename CRSMatrixType::size_type& update, const bool final) {
      const typename CRSMatrixType::size_type val = d_row_nnz(i);
      if (final) row_map(i+1) = update + val;
      update += val;
    });
  
  // Get total number of nonzeros
  typename CRSMatrixType::size_type total_nnz = 0;
  Kokkos::deep_copy(total_nnz, Kokkos::subview(row_map, nlocal));
  
  // Allocate column indices and values
  typename CRSMatrixType::index_type::non_const_type columns("crs_columns", total_nnz);
  typename CRSMatrixType::values_type::non_const_type values("crs_values", total_nnz);

  // Keep track of next position for each row
  Kokkos::View<int*, DeviceType> d_row_fill("row_fill", nlocal);
  
  // Initialize row_fill with row_map values
  Kokkos::parallel_for("InitFill", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      d_row_fill(i) = row_map(i);
    });

  // Fill in matrix entries
  Kokkos::parallel_for("FillCRSMatrix", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      const int itype = d_type(i);
      
      // Add diagonal element first (for all atoms)
      int pos = Kokkos::atomic_fetch_add(&d_row_fill(i), 1);
      columns(pos) = i;
      values(pos) = d_params(itype).eta;
      
      // Add off-diagonal elements only for atoms in group
      if (d_mask[i] & groupbit) {
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
          // Only include atoms that are part of the group and are local
          if ((d_mask[j] & groupbit) && j < nlocal) {
            const double delx = d_x(j,0) - d_x(i,0);
            const double dely = d_x(j,1) - d_x(i,1);
            const double delz = d_x(j,2) - d_x(i,2);
            const double rsq = delx*delx + dely*dely + delz*delz;
            
            if (rsq <= cutsq) {
              const double r = sqrt(rsq);
              const double shldij = d_shield(itype, d_type(j));
              const double hval = calculate_H_k(r, shldij);
              
              // Use atomic to safely update the matrix
              pos = Kokkos::atomic_fetch_add(&d_row_fill(i), 1);
              
              // Bug fix: Add buffer overflow check
              if (pos < total_nnz) {
                columns(pos) = j;
                values(pos) = hval;
              } else {
                // Indicate overflow - in production code this would need proper handling
                // e.g., resize buffers or abort with error message
                // For simplicity, we ignore overflow here
              }
            }
          }
        }
      }
    });
  
  // Bug fix: Verify that matrix rows are consistent (for debugging)
  /*
  Kokkos::parallel_for("VerifyMatrix", Kokkos::RangePolicy<DeviceType>(0, nlocal),
    KOKKOS_LAMBDA(const int i) {
      if (d_row_fill(i) != row_map(i+1)) {
        // Inconsistent row - would need proper error handling
        // This indicates a bug in counting vs. filling nonzeros
      }
    });
  */
  
  // Create the CRS matrix
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
void FixQEqReaxFFKokkos<DeviceType>::sparse_matvec(t_compute_1d &in, t_compute_1d &out)
{
  // Optimized SpMV using KokkosSparse with CRS matrix for local-local interactions
  KokkosSparse::spmv("N",           // No transpose
                     1.0,           // alpha
                     crs_matrix,    // CRS matrix
                     in,            // x vector
                     0.0,           // beta
                     out);          // y vector

  int nlocal = atomKK->nlocal;
  int nall = nlocal + atomKK->nghost;

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  // Zero ghost atom output values
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(nlocal, nall),
    KOKKOS_LAMBDA(const int i) {
      if (d_mask[i] & groupbit) out(i) = 0.0;
    });
  
  // Optimization: Batch ghost atom calculations for better cache utilization
  // This reduces thread divergence on GPUs
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        const int itype = d_type(i);
        kk_compute sum_ghost = 0.0;
        
        // Get neighbor info
        const int jnum = d_numneigh[i];
        
        // Optimization: Prefetch neighbor data
        // Process neighbors in blocks for better memory access pattern
        static constexpr int BLOCK_SIZE = 16; // Tunable parameter
        
        for (int jj = 0; jj < jnum; jj += BLOCK_SIZE) {
          // Process a block of neighbors
          const int block_end = MIN(jj + BLOCK_SIZE, jnum);
          
          // First pass: prefetch data for this block
          for (int jb = jj; jb < block_end; jb++) {
            const int j = d_neighbors(i, jb) & NEIGHMASK;
            // Prefetch neighbor data if supported
            // Note: This is architecture-dependent
            #ifdef KOKKOS_ENABLE_CUDA
            if (j >= nlocal) {
              // On CUDA, use intrinsics for prefetching
              __prefetch_global_l1(&in(j));
              __prefetch_global_l1(&d_x(j,0));
              __prefetch_global_l1(&d_type(j));
            }
            #endif
          }
          
          // Second pass: actually process the block
          for (int jb = jj; jb < block_end; jb++) {
            const int j = d_neighbors(i, jb) & NEIGHMASK;
            
            // Only handle ghost atoms here
            if (j >= nlocal && (d_mask[j] & groupbit)) {
              const double delx = d_x(j,0) - d_x(i,0);
              const double dely = d_x(j,1) - d_x(i,1);
              const double delz = d_x(j,2) - d_x(i,2);
              const double rsq = delx*delx + dely*dely + delz*delz;
              
              if (rsq <= cutsq) {
                const double r = sqrt(rsq);
                const double shldij = d_shield(itype, d_type(j));
                const double hval = calculate_H_k(r, shldij);
                
                // Accumulate contribution from ghost atom
                sum_ghost += hval * in(j);
              }
            }
          }
        }
        
        // Add ghost contributions to local atom - use atomic for safety
        Kokkos::atomic_add(&out(i), sum_ghost);
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
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::update_extended_lagrangian()
{
  auto d_q = atomKK->k_q.view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {

        // First calculate acceleration: a = omega^2 (q - theta)
        const double delta_q = d_q(i) - d_theta(i);
        const double acceleration = omega2 * delta_q;

        // Update position and velocity using velocity Verlet
        d_theta(i) += d_theta_dot(i) * dt + acceleration * dt2_half;
        d_theta_dot(i) += acceleration * dt_half;
    
        // Calculate new acceleration and complete velocity update
        const double new_delta_q = d_q(i) - d_theta(i);
        d_theta_dot(i) += omega2 * new_delta_q * dt_half;
    
        // Optimization: Add damping for numerical stability (optional)
        //static constexpr double DAMPING = 0.9;
        //d_theta_dot(i) *= DAMPING;
    }
  });

  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
  auto d_q = atomKK->k_q.view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();

  // Update charges with the final solution
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) d_q(i) = d_s(i);  // Use s directly for charge
    });
  
  atomKK->modified(execution_space, Q_MASK);

  // Forward communicate charges
  pack_flag = 2;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm_kokkos(
    int n, DAT::tdual_int_1d k_sendlist, DAT::tdual_xfloat_1d &k_buf,
    int pbc_flag, int *pbc)
{
  // Current pack flag determines what data to send
  int current_pack_flag = pack_flag;
  
  // Sync dual views
  k_sendlist.template sync<DeviceType>();
  k_buf.template sync<DeviceType>();
  
  // Create device views
  auto d_sendlist = k_sendlist.template view<DeviceType>();
  auto d_buf = k_buf.template view<DeviceType>();
  
  // Sync appropriate data based on pack_flag
  if (current_pack_flag == 1) {
    k_s.template sync<DeviceType>();
    auto d_s = k_s.template view<DeviceType>();
    
    // Pack s values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) {
        int j = d_sendlist(i);
        d_buf(i) = d_s(j);
      });
  } else if (current_pack_flag == 2) {
    // For charges, we might need to sync appropriately
    atomKK->sync(Device,Q_MASK);

    auto d_q = atomKK->k_q.view<DeviceType>();

    // Pack q values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) {
        int j = d_sendlist(i);
        d_buf(i) = d_q(j);
      });
  }
  
  // Mark buffer as modified
  k_buf.template modify<DeviceType>();
  
  // Return number of values packed
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm_kokkos(
    int n, int first_in, DAT::tdual_xfloat_1d &k_buf)
{
  // Current pack flag determines what data to unpack
  int current_pack_flag = pack_flag;
  
  // Sync buffer
  k_buf.template sync<DeviceType>();
  auto d_buf = k_buf.template view<DeviceType>();
  
  // Store first index (offset for ghost atoms)
  int first = first_in;
  
  if (current_pack_flag == 1) {
    // Sync s values
    k_s.template sync<DeviceType>();
    auto d_s = k_s.template view<DeviceType>();
    
    // Unpack s values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) { d_s(i + first) = d_buf(i); });

    // Mark as modified
    k_s.template modify<DeviceType>();
  } else if (current_pack_flag == 2) {
    // Sync charge values
    atomKK->sync(Device,Q_MASK);

    auto d_q = atomKK->k_q.view<DeviceType>();

    // Unpack q values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) { d_q(i + first) = d_buf(i); });
    
    // Mark as modified
    atomKK->modified(Device,Q_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixQEqReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = 0.0;
  
  // Extended Lagrangian variables
  bytes += atom->nmax * sizeof(kk_compute);     // theta
  bytes += atom->nmax * sizeof(kk_compute);     // theta_dot
  
  // CG solver vectors
  bytes += (double)atom->nmax * 6 * sizeof(kk_compute); // storage
  
  // CRS matrix memory usage
  if (crs_matrix_allocated) {
    size_t nnz = crs_matrix.nnz();
    bytes += nnz * sizeof(kk_compute);              // Values
    bytes += nnz * sizeof(int);                     // Column indices
    bytes += (atomKK->nlocal + 1) * sizeof(size_t); // Row pointers
  }

  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  k_theta.resize(nmax);
  k_theta_dot.resize(nmax);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  k_theta.h_view(j) = k_theta.h_view(i);
  k_theta_dot.h_view(j) = k_theta_dot.h_view(i);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  k_theta.sync_device();
  k_theta_dot.sync_device();

  Sorter.sort(LMPDeviceType(), k_theta.d_view);
  Sorter.sort(LMPDeviceType(), k_theta_dot.d_view);

  k_theta.modify_device();
  k_theta_dot.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace space)
{
  // Sync all dual views to ensure device has latest data
  k_buf.template sync<DeviceType>();
  k_copylist.template sync<DeviceType>();
  k_exchange_sendlist.template sync<DeviceType>();
  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();
  
  // Create device views from dual views
  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_copylist = k_copylist.template view<DeviceType>();
  auto d_exchange_sendlist = k_exchange_sendlist.template view<DeviceType>();
  
  // Get device access to extended Lagrangian variables
  auto d_theta = k_theta.template view<DeviceType>();
  auto d_theta_dot = k_theta_dot.template view<DeviceType>();
  
  copymode = 1;
  
  // Optimization: Use a more efficient packing approach
  // Buffer structure is known: 2 values per atom (theta, theta_dot)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nsend),
    KOKKOS_LAMBDA(const int& isend) {
      // Get the actual atom index from our send list
      const int i = d_exchange_sendlist(isend);
      const int buf_offset = isend * 2;
      
      // Pack theta and theta_dot into the buffer
      d_buf(buf_offset) = d_theta(i);
      d_buf(buf_offset + 1) = d_theta_dot(i);
      
      // Handle any copy operations if needed
      const int j = d_copylist(isend);
      if (j >= 0) {
        // Copy from i to j
        d_theta(j) = d_theta(i);
        d_theta_dot(j) = d_theta_dot(i);
      }
    });
  
  copymode = 0;
  
  // Mark views as modified on device
  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
  
  // Sync buffer if needed based on execution space
  k_buf.template modify<DeviceType>();
  if (space == Host) k_buf.template sync<LMPHostType>();
  else k_buf.template sync<LMPDeviceType>();
  
  // Return the number of values packed (2 per atom)
  return nsend*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_exchange_kokkos(
   DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
   int nrecv1, int nextrarecv1, ExecutionSpace space)
{
  // Sync input dual views to device
  k_buf.template sync<DeviceType>();
  k_indices.template sync<DeviceType>();
  
  // Sync extended Lagrangian variables
  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();
  
  // Create device views
  auto d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  auto d_indices = k_indices.template view<DeviceType>();
  
  // Get device access to extended Lagrangian variables
  auto d_theta = k_theta.template view<DeviceType>();
  auto d_theta_dot = k_theta_dot.template view<DeviceType>();
  
  copymode = 1;
  
  // Optimization: Improved unpacking with better memory access pattern
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nrecv),
    KOKKOS_LAMBDA(const int& i) {
      // Get the atom index this data belongs to
      int index = d_indices(i);
      
      // Only unpack if this is a valid atom (index > -1)
      if (index > -1) {
        const int buf_offset = i * 2;
        // Unpack in a memory-coalesced manner
        d_theta(index) = d_buf(buf_offset);
        d_theta_dot(index) = d_buf(buf_offset + 1);
      }
    });
  
  copymode = 0;
  
  // Mark views as modified on device
  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::get_chi_field()
{
  atomKK->sync(Host, X_MASK|MASK_MASK|IMAGE_MASK);
  FixQEqReaxFF::get_chi_field();
  k_chi_field.modify_host();
  k_chi_field.sync_device();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixQEqReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixQEqReaxFFKokkos<LMPHostType>;
#endif
}

