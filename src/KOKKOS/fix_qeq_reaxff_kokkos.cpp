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

static constexpr double EV_TO_KCAL_PER_MOL = 14.4;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxFFKokkos<DeviceType>::FixQEqReaxFFKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixQEqReaxFF(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | Q_MASK | MASK_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = X_MASK;

  nmax = 0;
  allocated_flag = 0;
  last_allocate = -1;

  // Extended Lagrangian parameters
  K = 2.0;  // recommended value from Niklasson papers
  omega = sqrt(K/(update->dt * update->dt));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxFFKokkos<DeviceType>::~FixQEqReaxFFKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_theta);
  memoryKK->destroy_kokkos(k_theta_dot);
  memoryKK->destroy_kokkos(k_chi_field, chi_field);
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

  last_allocate = -1;
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

  k_tap = DAT::tdual_float_1d("qeq/reaxff/kk:tap",8);
  d_tap = k_tap.template view<DeviceType>();

  for (i = 0; i < 8; i++)
    k_tap.h_view(i) = Tap[i];

  k_tap.template modify<LMPHostType>();
  k_tap.template sync<DeviceType>();
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

  d_q = atomKK->k_q.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  k_params.template sync<DeviceType>();
  k_shield.template sync<DeviceType>();
  k_tap.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  nn = list->inum;

  copymode = 1;

  // allocate
  allocate_array();

  // Build the sparse matrix directly in CRS format
  if (!allocated_flag || last_allocate < neighbor->lastcall) {
    build_crs_matrix();
    last_allocate = update->ntimestep;
  }

  // Initialize matvec
  if (efield) get_chi_field();

  // Initialize vectors for CG iteration
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqInitMatvec>(0, nn), *this);

  // Pack and communicate charges
  pack_flag = 1;
  k_s.template modify<DeviceType>();
  comm->forward_comm(this);
  k_s.template sync<DeviceType>();

  // Do one CG iteration using theta as initial guess
  one_cg_iter();

  // Calculate atomic charges
  calculate_q();

  // Update auxiliary variables using velocity Verlet
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqUpdate>(0, nn), *this);
  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();

  copymode = 0;

  if (!allocated_flag)
    allocated_flag = 1;

  atomKK->modified(execution_space, datamask_modify);

  d_neighbors = typename AT::t_neighbors_2d();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::allocate_array()
{

  utils::logmesg(lmp, "*** allocate_array() atom->nmax {} nmax {}\n", atom->nmax, nmax);

  if (atom->nmax > nmax) {
    nmax = atom->nmax;


    // Allocate regular QEq arrays
    memoryKK->create_kokkos(k_o,nmax,"qeq/reaxff/kk:k_o");
    d_o = k_o.template view<DeviceType>();

    d_r = t_compute_1d("qeq/reaxff/kk:r", nmax);
    d_p = t_compute_1d("qeq/reaxff/kk:p", nmax);
    d_d = t_compute_1d("qeq/reaxff/kk:d", nmax);
    d_Hdia_inv = t_compute_1d("qeq/reaxff/kk:Hdia_inv", nmax);
    d_b = t_compute_1d("qeq/reaxff/kk:b", nmax);

    k_s = tdual_compute_1d("qeq/reaxff/kk:s", nmax);
    d_s = k_s.template view<DeviceType>();

    memoryKK->create_kokkos(k_chi_field, chi_field, nmax, "qeq/reaxff/kk:chi_field");
    d_chi_field = k_chi_field.template view<DeviceType>();

    // Allocate extended Lagrangian variables
    memoryKK->create_kokkos(k_theta, nmax, "qeq/reaxff/kk:theta");
    memoryKK->create_kokkos(k_theta_dot, nmax, "qeq/reaxff/kk:theta_dot");
    d_theta = k_theta.template view<DeviceType>();
    d_theta_dot = k_theta_dot.template view<DeviceType>();

    for (int i = 0; i < nmax; i++) {
      k_theta.h_view(i) = 0.0;
      k_theta_dot.h_view(i) = 0.0;
    }

    k_theta.template modify<LMPHostType>();
    k_theta_dot.template modify<LMPHostType>();
  }

  // init_storage
  if (efield) get_chi_field();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxFFKokkos<DeviceType>::calculate_H_k(const double &r, const double &shld) const
{
  double taper, denom;

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
void FixQEqReaxFFKokkos<DeviceType>::build_crs_matrix()
{
  // Direct construction of CRS matrix for QEq
  int nlocal = atomKK->nlocal;

  // First, count interactions to allocate arrays
  Kokkos::View<int*, DeviceType> d_row_nnz("row_nnz", nlocal);

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  // Count nonzeros per row
  Kokkos::parallel_for("CountNNZ", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        int count = 1; // Diagonal element
        
        // Count off-diagonal elements
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
            if (rsq <= cutsq) count++;
          }
        }
        d_row_nnz(i) = count;
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

  // Fill in matrix entries using SparseRowView approach
  Kokkos::parallel_for("FillCRSMatrix", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        const int itype = d_type(i);
        
        // Calculate row start and end indices
        const int row_start = row_map(i);
        const int row_end = row_map(i+1);
        
        // Track current position in row
        int pos = row_start;
        
        // Add diagonal element first
        columns(pos) = i;
        values(pos) = d_params(itype).eta;
        pos++;
        
        // Add off-diagonal elements
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
              
              // Safety check to prevent buffer overrun
              if (pos < row_end) {
                columns(pos) = j;
                values(pos) = hval;
                pos++;
              }
            }
          }
        }
      }
    });
  
  // Create the CRS matrix - fixed constructor call with required 7 parameters
  crs_matrix = CRSMatrixType("H_crs", nlocal, nlocal, total_nnz, values, row_map, columns);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sparse_matvec(t_compute_1d &in, t_compute_1d &out)
{
  // Optimized SpMV using KokkosSparse with CRS matrix
  KokkosSparse::spmv("N",           // No transpose
                     1.0,           // alpha
                     crs_matrix,    // CRS matrix
                     in,            // x vector
                     0.0,           // beta
                     out);          // y vector

  int nlocal = atomKK->nlocal;

  // Handle ghost atoms separately - this is required because the CRS matrix only includes
  // local-local interactions, but we need local-ghost interactions as well
  int nall = nlocal + atomKK->nghost;

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  // Zero ghost atom output values
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(nlocal, nall),
    KOKKOS_LAMBDA(const int i) {
      if (d_mask[i] & groupbit) out(i) = 0.0;
    });
  
  // Handle interactions with ghost atoms
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        const int itype = d_type(i);
        kk_compute sum_ghost = 0.0;

        // Loop over neighbors
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
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
              
              // Add contribution from ghost atom
              sum_ghost += hval * in(j);
            }
          }
        }
        
        // Add ghost contributions to local atom
        out(i) += sum_ghost;
      }
    });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqInitMatvec, const int &ii) const
{
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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::one_cg_iter()
{
  // Simple function to perform one step of CG minimization

  // r = b - H*s
  sparse_matvec(d_s, d_o);  // H*s -> o
  
  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    comm->reverse_comm(this);
    k_o.template sync<DeviceType>();
  }

  // r = b - o (residual)
  // d = r * Hdia_inv (preconditioned residual)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (d_mask[i] & groupbit) {
        d_r(i) = d_b(i) - d_o(i);
        d_d(i) = d_r(i) * d_Hdia_inv[i];
      }
    });
  
  // Compute dot products using KokkosBlas
  double rnorm2 = KokkosBlas::dot(d_r, d_d);

  // Compute alpha = r*d / (d*H*d)
  sparse_matvec(d_d, d_o);  // H*d -> o

  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    comm->reverse_comm(this);
    k_o.template sync<DeviceType>();
  }

  // Compute dHd = d·o
  double dHd = KokkosBlas::dot(d_d, d_o);

  // alpha = rnorm2 / dHd
  double alpha = rnorm2 / dHd;

  // s = s + alpha*d
  KokkosBlas::axpy(alpha, d_d, d_s);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqUpdate, const int &ii) const
{
  const int i = d_ilist[ii];
  if (d_mask[i] & groupbit) {
    // Velocity Verlet integration for auxiliary variables
    const double dt = update->dt;
    
    // First calculate acceleration: a = omega^2 (q - theta)
    const double acc = omega * omega * (d_q(i) - d_theta(i));

    // Update position and velocity
    d_theta(i) += d_theta_dot(i) * dt + 0.5 * acc * dt * dt;
    d_theta_dot(i) += 0.5 * acc * dt;
    
    // Store the second half of velocity update for the next step
    // Calculate new acceleration based on updated position
    const double new_acc = omega * omega * (d_q(i) - d_theta(i));
    d_theta_dot(i) += 0.5 * new_acc * dt;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
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

    // Unpack q values
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, n),
      KOKKOS_LAMBDA(const int& i) { d_q(i + first) = d_buf(i); });
    
    // Mark as modified
    atomKK->modified(Device,Q_MASK);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixQEqReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = atom->nmax * sizeof(kk_compute);     // theta
  bytes += atom->nmax * sizeof(kk_compute);           // theta_dot
  bytes += (double)atom->nmax*6 * sizeof(kk_compute); // storage
  
  // CRS matrix memory usage
  size_t nnz = crs_matrix.nnz();
  bytes += nnz * sizeof(kk_compute);              // Values
  bytes += nnz * sizeof(int);                     // Column indices
  bytes += (atomKK->nlocal + 1) * sizeof(size_t); // Row pointers

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

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

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

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
  
  // Since we know exactly where each atom's data goes (fixed size of 2 per atom),
  // a simple parallel_for is sufficient
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nsend),
    KOKKOS_LAMBDA(const int& mysend) {
      // Get the actual atom index from our send list
      const int i = d_exchange_sendlist(mysend);
      
      // Pack theta and theta_dot into the buffer
      d_buf(mysend*2) = d_theta(i);
      d_buf(mysend*2+1) = d_theta_dot(i);
      
      // Handle any copy operations if needed
      const int j = d_copylist(mysend);
      // If j > -1, this was potentially a copy operation in the original code
      // However, no copy logic was implemented in the original
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
  
  // Process received data
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nrecv),
    KOKKOS_LAMBDA(const int& i) {
      // Get the atom index this data belongs to
      int index = d_indices(i);
      
      // Only unpack if this is a valid atom (index > -1)
      if (index > -1) {
        // For QEq, we have 2 values per atom at fixed positions
        // Each atom's data starts at position i*2
        d_theta(index) = d_buf(i*2);
        d_theta_dot(index) = d_buf(i*2+1);
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
