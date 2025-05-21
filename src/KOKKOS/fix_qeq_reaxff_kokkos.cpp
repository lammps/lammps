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

#include "fix_qeq_reaxff_schur_kokkos.hpp"

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

#include "fix_qeq_reaxff_comm_kokkos.hpp"

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

namespace LAMMPS_NS {
template class FixQEqReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixQEqReaxFFKokkos<LMPHostType>;
#endif
}

