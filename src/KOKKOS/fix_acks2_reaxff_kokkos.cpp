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
   Contributing authors: Stan Moore (SNL)
                         Mitch Murphy (alphataubio at gmail)
                         - add target_charge option
                         - sparse BSR mixed precision Schur CG
------------------------------------------------------------------------- */

#include "fix_acks2_reaxff_kokkos.h"

#include "atom.h"
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
FixACKS2ReaxFFKokkos<DeviceType>::
FixACKS2ReaxFFKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixACKS2ReaxFF(lmp, narg, arg)
{
  kokkosable = 1;
  sort_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | MASK_MASK | Q_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = Q_MASK;

  nmax = 0;
  allocated_flag = 0;

  d_mfill_offset = typename AT::t_bigint_scalar("acks2/kk:mfill_offset");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixACKS2ReaxFFKokkos<DeviceType>::~FixACKS2ReaxFFKokkos()
{
  if (copymode) return;
  deallocate_array();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::post_constructor()
{
  grow_arrays(atom->nmax);
  pertype_parameters(pertype_option);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init()
{
  atomKK->k_q.modify<LMPHostType>();
  atomKK->k_q.sync<DeviceType>();

  FixACKS2ReaxFF::init();

  // Set up neighbor list request for KOKKOS
  neighflag = lmp->kokkos->neighflag_qeq;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  // Set up parameter arrays
  int ntypes = atom->ntypes;
  k_params = Kokkos::DualView<params_acks2*,Kokkos::LayoutRight,DeviceType>
    ("FixACKS2ReaxFF::params",ntypes+1);
  params = k_params.template view<DeviceType>();

  for (int n = 1; n <= ntypes; n++) {
    k_params.h_view(n).chi = chi[n];
    k_params.h_view(n).eta = eta[n];
    k_params.h_view(n).gamma = gamma[n];
    k_params.h_view(n).bcut_acks2 = bcut_acks2[n];
  }
  k_params.template modify<LMPHostType>();

  cutsq = swb * swb;

  init_shielding_k();
  init_bsr();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init_shielding_k()
{
  int i,j;
  int ntypes = atom->ntypes;

  k_shield = DAT::tdual_float_2d("acks2/kk:shield",ntypes+1,ntypes+1);
  d_shield = k_shield.template view<DeviceType>();
  k_shield.template modify<DeviceType>();

  k_bcut = DAT::tdual_float_2d("acks2/kk:bcut",ntypes+1,ntypes+1);
  d_bcut = k_bcut.template view<DeviceType>();
  k_bcut.template modify<DeviceType>();

  k_tap = DAT::tdual_float_1d("acks2/kk:tap",8);
  d_tap = k_tap.template view<DeviceType>();

  for (i = 1; i <= ntypes; ++i) {
    for (j = 1; j <= ntypes; ++j) {
      k_shield.h_view(i,j) = shld[i][j];
      k_bcut.h_view(i,j) = bcut[i][j];
    }
  }

  for (i = 0; i < 8; i++) {
    k_tap.h_view(i) = Tap[i];
  }

  k_shield.template modify<LMPHostType>();
  k_shield.template sync<DeviceType>();

  k_bcut.template modify<LMPHostType>();
  k_bcut.template sync<DeviceType>();

  k_tap.template modify<LMPHostType>();
  k_tap.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init_bsr()
{
  // No history initialization needed for BSR + Schur CG
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::allocate_array()
{
  if (allocated_flag == 0) {
    nmax = atom->nmax;
    
    k_s = DAT::tdual_float_1d("acks2:s",2*nmax+2);
    d_s = k_s.view<DeviceType>();

    k_X_diag = DAT::tdual_float_1d("acks2:X_diag",nmax);
    d_X_diag = k_X_diag.view<DeviceType>();

    k_chi_field = DAT::tdual_float_1d("acks2:chi_field",nmax);
    d_chi_field = k_chi_field.view<DeviceType>();

    k_b_s = DAT::tdual_float_1d("acks2:b_s",2*nmax+2);
    d_b_s = k_b_s.view<DeviceType>();

    allocated_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::deallocate_array()
{
  if (allocated_flag == 0) return;
  allocated_flag = 0;
}

/* ---------------------------------------------------------------------- */

#include "fix_acks2_reaxff_schur_kokkos.hpp"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixACKS2ReaxFFKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                         int /* pbc_flag */, int * /* pbc */)
{
  int m = 0;

  if (pack_flag == 1) {
    k_chi_field.sync_host();
    for (m = 0; m < n; m++) {
      buf[m] = k_chi_field.h_view[list[m]];
    }
  } else if (pack_flag == 2) {
    k_s.sync_host();
    for (m = 0; m < n; m++) {
      buf[m] = k_s.h_view[list[m]];
      buf[n+m] = k_s.h_view[list[m]+nmax];
    }
    m = 2;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1) {
    k_chi_field.sync_host();
    for (m = 0, i = first; m < n; m++, i++) {
      k_chi_field.h_view[i] = buf[m];
    }
    k_chi_field.modify_host();
  } else if (pack_flag == 2) {
    k_s.sync_host();
    for (m = 0, i = first; m < n; m++, i++) {
      k_s.h_view[i] = buf[m];
      k_s.h_view[i+nmax] = buf[n+m];
    }
    k_s.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixACKS2ReaxFFKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  if (pack_flag == 4) {
    k_X_diag.sync_host();
    for (m = 0, i = first; m < n; m++, i++) {
      buf[m] = k_X_diag.h_view[i];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m;

  if (pack_flag == 4) {
    k_X_diag.sync_host();
    for (m = 0; m < n; m++) {
      k_X_diag.h_view[list[m]] += buf[m];
    }
    k_X_diag.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixACKS2ReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = 0.0;

  int n = atom->nmax;
  int N = nlocal + atom->nghost;

  // Basic arrays
  bytes += (double)n * sizeof(double);        // X_diag
  bytes += (double)n * sizeof(double);        // chi_field
  bytes += (double)2 * n * sizeof(double);    // b_s, s (solution vector)

  // BSR matrix storage
  if (allocated_flag) {
    // Get BSR matrix sizes
    int nnz_blocks = bsr_matrix.graph.entries.extent(0);
    int num_block_rows = bsr_matrix.numRows();
    
    bytes += (double)nnz_blocks * 4 * sizeof(double);     // values (4 per block)
    bytes += (double)nnz_blocks * sizeof(int);            // block_cols
    bytes += (double)(num_block_rows + 1) * sizeof(int);  // block_row_map
    bytes += (double)num_block_rows * sizeof(int);        // row_block_counts
  }

  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::get_chi_field()
{
  atomKK->sync(Host,X_MASK|MASK_MASK|IMAGE_MASK);
  FixACKS2ReaxFF::get_chi_field();
  k_chi_field.modify_host();
  k_chi_field.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::allocate_matrix()
{
  // BSR matrix allocation is handled by allocate_bsr_matrix()
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::deallocate_matrix()
{
  // BSR matrix deallocation happens automatically when views go out of scope
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::setup_pre_force(int /* vflag */)
{
  neighbor->build_one(list);
  
  // Ensure proper allocation
  allocate_array();
  if (!allocated_flag) {
    allocate_bsr_matrix();
    allocated_flag = 1;
  }
  
  // Build matrix structure
  build_bsr_matrix();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixACKS2ReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixACKS2ReaxFFKokkos<LMPHostType>;
#endif
}
