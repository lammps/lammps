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
     - Implement extended lagrangian (Nomura, 2015)
       https://doi.org/10.1016/j.cpc.2015.02.023
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

  memoryKK->destroy_kokkos(k_theta, nullptr);
  memoryKK->destroy_kokkos(k_theta_dot, nullptr);
  memoryKK->destroy_kokkos(k_chi_field, chi_field);
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
  k_params = Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType>
    ("FixQEqReaxFF::params", ntypes+1);
  params = k_params.template view<DeviceType>();

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

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  nlocal = atomKK->nlocal;
  newton_pair = force->newton_pair;

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

  // Build the sparse matrix directly in CSR format
  if (!allocated_flag || last_allocate < neighbor->lastcall) {
    build_csr_matrix();
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
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::num_neigh_item(int ii, bigint &totneigh) const
{
  const int i = d_ilist[ii];
  totneigh += d_numneigh[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::allocate_array()
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;

    // Allocate regular QEq arrays
    k_o = DAT::tdual_float_1d("qeq/reaxff/kk:o", nmax);
    d_o = t_compute_1d(k_o.template view<DeviceType>().data(), nmax);
    h_o = k_o.h_view;

    d_r = t_compute_1d("qeq/reaxff/kk:r", nmax);
    d_p = t_compute_1d("qeq/reaxff/kk:p", nmax);
    d_d = t_compute_1d("qeq/reaxff/kk:d", nmax);
    d_Hdia_inv = t_compute_1d("qeq/reaxff/kk:Hdia_inv", nmax);
    d_b = t_compute_1d("qeq/reaxff/kk:b", nmax);

    k_s = DAT::tdual_float_1d("qeq/reaxff/kk:s", nmax);
    d_s = t_compute_1d(k_s.template view<DeviceType>().data(), nmax);
    h_s = k_s.h_view;

    memoryKK->create_kokkos(k_chi_field, chi_field, nmax, "qeq/reaxff/kk:chi_field");
    d_chi_field = k_chi_field.template view<DeviceType>();

    // Allocate extended Lagrangian variables
    memoryKK->create_kokkos(k_theta, nullptr, nmax, "qeq/reaxff/kk:theta");
    d_theta = t_compute_1d(k_theta.template view<DeviceType>().data(), nmax);

    memoryKK->create_kokkos(k_theta_dot, nullptr, nmax, "qeq/reaxff/kk:theta_dot");
    d_theta_dot = t_compute_1d(k_theta_dot.template view<DeviceType>().data(), nmax);

    // Initialize auxiliary variables
    k_theta.template sync<LMPHostType>();
    k_theta_dot.template sync<LMPHostType>();
    
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
void FixQEqReaxFFKokkos<DeviceType>::build_csr_matrix()
{
  // Direct construction of CSR matrix for QEq
  int num_rows = nlocal;
  
  // First, count interactions to allocate arrays
  Kokkos::View<int*, DeviceType> d_row_nnz("row_nnz", num_rows);
  
  // Count nonzeros per row
  Kokkos::parallel_for("CountNNZ", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (mask[i] & groupbit) {
        int count = 1; // Diagonal element
        
        // Count off-diagonal elements
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
          // Only include atoms that are part of the group and are local
          if ((mask[j] & groupbit) && j < nlocal) {
            const double delx = x(j,0) - x(i,0);
            const double dely = x(j,1) - x(i,1);
            const double delz = x(j,2) - x(i,2);
            const double rsq = delx*delx + dely*dely + delz*delz;
            
            if (rsq <= cutsq) {
              count++;
            }
          }
        }
        d_row_nnz(i) = count;
      }
    });
  
  // Create row map from row_nnz (exclusive scan)
  typename CSRMatrixType::row_map_type::non_const_type row_map("csr_row_map", num_rows + 1);
  
  // Set first element to 0
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, 1),
    KOKKOS_LAMBDA(const int) {
      row_map(0) = 0;
    });
  
  // Exclusive scan to compute row offsets
  Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType>(0, num_rows),
    KOKKOS_LAMBDA(const int i, typename CSRMatrixType::size_type& update, const bool final) {
      const typename CSRMatrixType::size_type val = d_row_nnz(i);
      if (final) {
        row_map(i+1) = update + val;
      }
      update += val;
    });
  
  // Get total number of nonzeros
  typename CSRMatrixType::size_type total_nnz = 0;
  Kokkos::deep_copy(total_nnz, Kokkos::subview(row_map, num_rows));
  
  // Allocate column indices and values
  typename CSRMatrixType::index_type::non_const_type columns("csr_columns", total_nnz);
  typename CSRMatrixType::values_type::non_const_type values("csr_values", total_nnz);
  
  // Fill in matrix entries
  Kokkos::View<int*, DeviceType> d_nnz_count("nnz_count", num_rows);
  Kokkos::deep_copy(d_nnz_count, 0);
  
  Kokkos::parallel_for("FillCSRMatrix", Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (mask[i] & groupbit) {
        const int itype = type(i);
        
        // Add diagonal element first
        const int local_col = Kokkos::atomic_fetch_add(&d_nnz_count(i), 1);
        const int pos = row_map(i) + local_col;
        columns(pos) = i;
        values(pos) = params(itype).eta;
        
        // Add off-diagonal elements
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
          // Only include atoms that are part of the group and are local
          if ((mask[j] & groupbit) && j < nlocal) {
            const double delx = x(j,0) - x(i,0);
            const double dely = x(j,1) - x(i,1);
            const double delz = x(j,2) - x(i,2);
            const double rsq = delx*delx + dely*dely + delz*delz;
            
            if (rsq <= cutsq) {
              const double r = sqrt(rsq);
              const double shldij = d_shield(itype, type(j));
              const double hval = calculate_H_k(r, shldij);
              
              const int local_col = Kokkos::atomic_fetch_add(&d_nnz_count(i), 1);
              const int pos = row_map(i) + local_col;
              columns(pos) = j;
              values(pos) = hval;
            }
          }
        }
      }
    });
  
  // Create the CSR matrix
  csr_matrix = CSRMatrixType("H_csr", num_rows, num_rows, values, row_map, columns);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sparse_matvec(t_compute_1d &in, t_compute_1d &out)
{
  // Optimized SpMV using KokkosSparse with CSR matrix
  KokkosSparse::spmv("N",           // No transpose
                     1.0,           // alpha
                     csr_matrix,    // CSR matrix
                     in,            // x vector
                     0.0,           // beta
                     out);          // y vector
  
  // Handle ghost atoms separately - this is required because the CSR matrix only includes
  // local-local interactions, but we need local-ghost interactions as well
  int nall = atom->nlocal + atom->nghost;
  
  // Zero ghost atom output values
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(atom->nlocal, nall),
    KOKKOS_LAMBDA(const int i) {
      if (mask[i] & groupbit) {
        out(i) = 0.0;
      }
    });
  
  // Handle interactions with ghost atoms
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (mask[i] & groupbit) {
        const int itype = type(i);
        compute_float sum_ghost = 0.0;
        
        // Loop over neighbors
        const int jnum = d_numneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
          int j = d_neighbors(i, jj);
          j &= NEIGHMASK;
          
          // Only handle ghost atoms here
          if (j >= nlocal && (mask[j] & groupbit)) {
            const double delx = x(j,0) - x(i,0);
            const double dely = x(j,1) - x(i,1);
            const double delz = x(j,2) - x(i,2);
            const double rsq = delx*delx + dely*dely + delz*delz;
            
            if (rsq <= cutsq) {
              const double r = sqrt(rsq);
              const double shldij = d_shield(itype, type(j));
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
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b[i] = -params(itype).chi - d_chi_field[i];
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
      if (mask[i] & groupbit) {
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
  if (mask[i] & groupbit) {
    // Velocity Verlet integration for auxiliary variables
    const double dt = update->dt;
    
    // First calculate acceleration: a = omega^2 (q - theta)
    const double acc = omega * omega * (q(i) - d_theta(i));
    
    // Update position and velocity
    d_theta(i) += d_theta_dot(i) * dt + 0.5 * acc * dt * dt;
    d_theta_dot(i) += 0.5 * acc * dt;
    
    // Store the second half of velocity update for the next step
    // Calculate new acceleration based on updated position
    const double new_acc = omega * omega * (q(i) - d_theta(i));
    d_theta_dot(i) += 0.5 * new_acc * dt;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
  // Update charges with the final solution
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqCalculateQ>(0, nn),
    KOKKOS_LAMBDA(const int ii) {
      const int i = d_ilist[ii];
      if (mask[i] & groupbit) {
        q(i) = d_s(i);  // Use s directly for charge
      }
    });
  
  atomKK->modified(execution_space, Q_MASK);

  // Forward communicate charges
  pack_flag = 2;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqCalculateQ, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) q(i) = d_s(i);  // Use s directly for charge
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                        DAT::tdual_xfloat_1d &k_buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  d_buf = k_buf.view<DeviceType>();
  
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqPackForwardComm>(0, n), *this);
    
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqPackForwardComm, const int &i) const {
  int j = d_sendlist(i);

  if (pack_flag == 1) d_buf[i] = d_s(j);
  else if (pack_flag == 2) d_buf[i] = q[j];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  d_buf = buf.view<DeviceType>();
  
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqUnpackForwardComm>(0, n), *this);

  if (pack_flag == 2)
    atomKK->modified(execution_space, Q_MASK); // needed for auto_sync
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqUnpackForwardComm, const int &i) const {
  if (pack_flag == 1) d_s(i+first) = d_buf[i];
  else if (pack_flag == 2) q[i + first] = d_buf[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  if (pack_flag == 1) {
    k_s.template sync<LMPHostType>();
    for (m = 0; m < n; m++) buf[m] = h_s(list[m]);
  } else if (pack_flag == 2) {
    atomKK->sync(Host, Q_MASK);
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  }

  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1) {
    k_s.template sync<LMPHostType>();
    for (m = 0, i = first; m < n; m++, i++) h_s(i) = buf[m];
    k_s.template modify<LMPHostType>();
  } else if (pack_flag == 2) {
    atomKK->sync(Host, Q_MASK);
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
    atomKK->modified(Host, Q_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  k_o.sync_host();
  for (m = 0, i = first; m < n; m++, i++) {
    buf[m] = h_o(i);
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_o.sync_host();
  for (int m = 0; m < n; m++) {
    h_o(list[m]) += buf[m];
  }
  k_o.modify_host();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixQEqReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = atom->nmax * sizeof(compute_float); // theta
  bytes += atom->nmax * sizeof(compute_float);       // theta_dot
  bytes += (double)atom->nmax*6 * sizeof(compute_float); // storage
  
  // CSR matrix memory usage
  size_t nnz = csr_matrix.nnz();
  bytes += nnz * sizeof(compute_float);              // Values
  bytes += nnz * sizeof(int);                        // Column indices
  bytes += (nlocal + 1) * sizeof(size_t);            // Row pointers

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::grow_arrays(int nmax)
{
  // Allocate auxiliary variables
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();
  
  k_theta.template modify<LMPHostType>(); // force reallocation on host
  k_theta_dot.template modify<LMPHostType>();
  
  memoryKK->grow_kokkos(k_theta, nullptr, nmax, "qeq/reaxff/kk:theta");
  memoryKK->grow_kokkos(k_theta_dot, nullptr, nmax, "qeq/reaxff/kk:theta_dot");
  
  d_theta = t_compute_1d(k_theta.template view<DeviceType>().data(), nmax);
  d_theta_dot = t_compute_1d(k_theta_dot.template view<DeviceType>().data(), nmax);
  
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
void FixQEqReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyValueType, BinOp> &Sorter)
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
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqPackExchange, const int &mysend) const {
  const int i = d_exchange_sendlist(mysend);

  d_buf(mysend*2) = d_theta(i);
  d_buf(mysend*2+1) = d_theta_dot(i);

  const int j = d_copylist(mysend);

  if (j > -1) {
    d_theta(i) = d_theta(j);
    d_theta_dot(i) = d_theta_dot(j);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_exchange_kokkos(
   const int &nsend, DAT::tdual_xfloat_2d &k_buf,
   DAT::tdual_int_1d k_exchange_sendlist, DAT::tdual_int_1d k_copylist,
   ExecutionSpace /*space*/)
{
  k_buf.sync<DeviceType>();
  k_copylist.sync<DeviceType>();
  k_exchange_sendlist.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_copylist = k_copylist.view<DeviceType>();
  d_exchange_sendlist = k_exchange_sendlist.view<DeviceType>();
  this->nsend = nsend;

  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqPackExchange>(0,nsend),*this);

  copymode = 0;

  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();

  return nsend*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqUnpackExchange, const int &i) const
{
  int index = d_indices(i);

  if (index > -1) {
    d_theta(index) = d_buf(i*2);
    d_theta_dot(index) = d_buf(i*2+1);
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf, DAT::tdual_int_1d &k_indices, int nrecv,
  int /*nrecv1*/, int /*nextrarecv1*/,
  ExecutionSpace /*space*/)
{
  k_buf.sync<DeviceType>();
  k_indices.sync<DeviceType>();

  d_buf = typename ArrayTypes<DeviceType>::t_xfloat_1d_um(
    k_buf.template view<DeviceType>().data(),
    k_buf.extent(0)*k_buf.extent(1));
  d_indices = k_indices.view<DeviceType>();

  k_theta.template sync<DeviceType>();
  k_theta_dot.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqUnpackExchange>(0,nrecv),*this);

  copymode = 0;

  k_theta.template modify<DeviceType>();
  k_theta_dot.template modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  buf[0] = k_theta.h_view(i);
  buf[1] = k_theta_dot.h_view(i);

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();

  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_theta.template sync<LMPHostType>();
  k_theta_dot.template sync<LMPHostType>();

  k_theta.h_view(nlocal) = buf[0];
  k_theta_dot.h_view(nlocal) = buf[1];

  k_theta.template modify<LMPHostType>();
  k_theta_dot.template modify<LMPHostType>();

  return 2;
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
