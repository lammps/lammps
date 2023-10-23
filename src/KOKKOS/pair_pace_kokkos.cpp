// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   aE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_pace_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace-evaluator/ace_radial.h"
#include <cstring>

namespace LAMMPS_NS {
struct ACEImpl {
  ACEImpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACEImpl()
  {
    delete basis_set;
    delete ace;
  }
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

enum{FS,FS_SHIFTEDSCALED};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairPACEKokkos<DeviceType>::PairPACEKokkos(LAMMPS *lmp) : PairPACE(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairPACEKokkos<DeviceType>::~PairPACEKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  // deallocate views of views in serial to prevent issues in Kokkos tools

  if (k_splines_gk.h_view.data()) {
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        k_splines_gk.h_view(i, j).deallocate();
        k_splines_rnl.h_view(i, j).deallocate();
        k_splines_hc.h_view(i, j).deallocate();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::grow(int natom, int maxneigh)
{
  auto basis_set = aceimpl->basis_set;

  if ((int)A.extent(0) < natom) {

    MemKK::realloc_kokkos(A, "pace:A", natom, nelements, nradmax + 1, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(A_rank1, "pace:A_rank1", natom, nelements, nradbase);

    MemKK::realloc_kokkos(A_list, "pace:A_list", natom, idx_rho_max, basis_set->rankmax);
    //size is +1 of max to avoid out-of-boundary array access in double-triangular scheme
    MemKK::realloc_kokkos(A_forward_prod, "pace:A_forward_prod", natom, idx_rho_max, basis_set->rankmax + 1);

    MemKK::realloc_kokkos(e_atom, "pace:e_atom", natom);
    MemKK::realloc_kokkos(rhos, "pace:rhos", natom, basis_set->ndensitymax + 1); // +1 density for core repulsion
    MemKK::realloc_kokkos(dF_drho, "pace:dF_drho", natom, basis_set->ndensitymax + 1); // +1 density for core repulsion

    MemKK::realloc_kokkos(weights, "pace:weights", natom, nelements, nradmax + 1, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(weights_rank1, "pace:weights_rank1", natom, nelements, nradbase);

    // hard-core repulsion
    MemKK::realloc_kokkos(rho_core, "pace:rho_core", natom);
    MemKK::realloc_kokkos(dF_drho_core, "pace:dF_drho_core", natom);
    MemKK::realloc_kokkos(dB_flatten, "pace:dB_flatten", natom, idx_rho_max, basis_set->rankmax);
  }

  if (((int)ylm.extent(0) < natom) || ((int)ylm.extent(1) < maxneigh)) {

    // radial functions
    MemKK::realloc_kokkos(fr, "pace:fr", natom, maxneigh, nradmax, lmax + 1);
    MemKK::realloc_kokkos(dfr, "pace:dfr", natom, maxneigh, nradmax, lmax + 1);
    MemKK::realloc_kokkos(gr, "pace:gr", natom, maxneigh, nradbase);
    MemKK::realloc_kokkos(dgr, "pace:dgr", natom, maxneigh, nradbase);
    const int max_num_functions = MAX(nradbase, nradmax*(lmax + 1));
    MemKK::realloc_kokkos(d_values, "pace:d_values", natom, maxneigh, max_num_functions);
    MemKK::realloc_kokkos(d_derivatives, "pace:d_derivatives", natom, maxneigh, max_num_functions);

    // hard-core repulsion
    MemKK::realloc_kokkos(cr, "pace:cr", natom, maxneigh);
    MemKK::realloc_kokkos(dcr, "pace:dcr", natom, maxneigh);

    // spherical harmonics
    MemKK::realloc_kokkos(plm, "pace:plm", natom, maxneigh, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(dplm, "pace:dplm", natom, maxneigh, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(ylm, "pace:ylm", natom, maxneigh, (lmax + 1) * (lmax + 1));
    MemKK::realloc_kokkos(dylm, "pace:dylm", natom, maxneigh, (lmax + 1) * (lmax + 1));

    // short neigh list
    MemKK::realloc_kokkos(d_ncount, "pace:ncount", natom);
    MemKK::realloc_kokkos(d_mu, "pace:mu", natom, maxneigh);
    MemKK::realloc_kokkos(d_rhats, "pace:rhats", natom, maxneigh);
    MemKK::realloc_kokkos(d_rnorms, "pace:rnorms", natom, maxneigh);
    MemKK::realloc_kokkos(d_nearest, "pace:nearest", natom, maxneigh);

    MemKK::realloc_kokkos(f_ij, "pace:f_ij", natom, maxneigh);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_pertype()
{
  auto basis_set = aceimpl->basis_set;

  MemKK::realloc_kokkos(d_rho_core_cutoff, "pace:rho_core_cutoff", nelements);
  MemKK::realloc_kokkos(d_drho_core_cutoff, "pace:drho_core_cutoff", nelements);
  MemKK::realloc_kokkos(d_E0vals, "pace:E0vals", nelements);
  MemKK::realloc_kokkos(d_ndensity, "pace:ndensity", nelements);
  MemKK::realloc_kokkos(d_npoti, "pace:npoti", nelements);

  auto h_rho_core_cutoff = Kokkos::create_mirror_view(d_rho_core_cutoff);
  auto h_drho_core_cutoff = Kokkos::create_mirror_view(d_drho_core_cutoff);
  auto h_E0vals = Kokkos::create_mirror_view(d_E0vals);
  auto h_ndensity = Kokkos::create_mirror_view(d_ndensity);
  auto h_npoti = Kokkos::create_mirror_view(d_npoti);

  for (int n = 0; n < nelements; n++) {
    h_rho_core_cutoff[n] = basis_set->map_embedding_specifications.at(n).rho_core_cutoff;
    h_drho_core_cutoff[n] = basis_set->map_embedding_specifications.at(n).drho_core_cutoff;

    h_E0vals(n)= basis_set->E0vals(n);

    h_ndensity(n) = basis_set->map_embedding_specifications.at(n).ndensity;

    string npoti = basis_set->map_embedding_specifications.at(n).npoti;
    if (npoti == "FinnisSinclair")
      h_npoti(n) = FS;
    else if (npoti == "FinnisSinclairShiftedScaled")
      h_npoti(n) = FS_SHIFTEDSCALED;
  }

  Kokkos::deep_copy(d_rho_core_cutoff, h_rho_core_cutoff);
  Kokkos::deep_copy(d_drho_core_cutoff, h_drho_core_cutoff);
  Kokkos::deep_copy(d_E0vals, h_E0vals);
  Kokkos::deep_copy(d_ndensity, h_ndensity);
  Kokkos::deep_copy(d_npoti, h_npoti);

  MemKK::realloc_kokkos(d_wpre, "pace:wpre", nelements, basis_set->ndensitymax);
  MemKK::realloc_kokkos(d_mexp, "pace:mexp", nelements, basis_set->ndensitymax);

  auto h_wpre = Kokkos::create_mirror_view(d_wpre);
  auto h_mexp = Kokkos::create_mirror_view(d_mexp);

  for (int n = 0; n < nelements; n++) {
    const int ndensity = basis_set->map_embedding_specifications.at(n).ndensity;
    for (int p = 0; p < ndensity; p++) {
      h_wpre(n, p) = basis_set->map_embedding_specifications.at(n).FS_parameters.at(p * 2 + 0);
      h_mexp(n, p) = basis_set->map_embedding_specifications.at(n).FS_parameters.at(p * 2 + 1);
    }
  }

  Kokkos::deep_copy(d_wpre, h_wpre);
  Kokkos::deep_copy(d_mexp, h_mexp);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_splines()
{
  auto basis_set = aceimpl->basis_set;

  if (k_splines_gk.d_view.data()) {
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        k_splines_gk.h_view(i, j).deallocate();
        k_splines_rnl.h_view(i, j).deallocate();
        k_splines_hc.h_view(i, j).deallocate();
      }
    }
  }

  k_splines_gk = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_gk", nelements, nelements);
  k_splines_rnl = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_rnl", nelements, nelements);
  k_splines_hc = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_hc", nelements, nelements);

  ACERadialFunctions* radial_functions = dynamic_cast<ACERadialFunctions*>(basis_set->radial_functions);

  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      k_splines_gk.h_view(i, j) = radial_functions->splines_gk(i, j);
      k_splines_rnl.h_view(i, j) = radial_functions->splines_rnl(i, j);
      k_splines_hc.h_view(i, j) = radial_functions->splines_hc(i, j);
    }
  }

  k_splines_gk.modify_host();
  k_splines_rnl.modify_host();
  k_splines_hc.modify_host();

  k_splines_gk.sync_device();
  k_splines_rnl.sync_device();
  k_splines_hc.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_tilde()
{
  auto basis_set = aceimpl->basis_set;

  // flatten loops, get per-element count and max

  idx_rho_max = 0;
  int total_basis_size_max = 0;

  MemKK::realloc_kokkos(d_idx_rho_count, "pace:idx_rho_count", nelements);
  auto h_idx_rho_count = Kokkos::create_mirror_view(d_idx_rho_count);

  for (int n = 0; n < nelements; n++) {
    int idx_rho = 0;
    const int total_basis_size_rank1 = basis_set->total_basis_size_rank1[n];
    const int total_basis_size = basis_set->total_basis_size[n];

    ACECTildeBasisFunction *basis = basis_set->basis[n];

    // rank=1
    for (int func_rank1_ind = 0; func_rank1_ind < total_basis_size_rank1; ++func_rank1_ind)
      idx_rho++;

    // rank > 1
    for (int func_ind = 0; func_ind < total_basis_size; ++func_ind) {
      ACECTildeBasisFunction *func = &basis[func_ind];

      // loop over {ms} combinations in sum
      for (int ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind)
        idx_rho++;
    }
    h_idx_rho_count(n) = idx_rho;
    idx_rho_max = MAX(idx_rho_max, idx_rho);
    total_basis_size_max = MAX(total_basis_size_max, total_basis_size_rank1 + total_basis_size);
  }

  Kokkos::deep_copy(d_idx_rho_count, h_idx_rho_count);

  MemKK::realloc_kokkos(d_rank, "pace:rank", nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_num_ms_combs, "pace:num_ms_combs", nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_offsets, "pace:offsets", nelements, idx_rho_max);
  MemKK::realloc_kokkos(d_mus, "pace:mus", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ns, "pace:ns", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ls, "pace:ls", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ms_combs, "pace:ms_combs", nelements, idx_rho_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ctildes, "pace:ctildes", nelements, idx_rho_max, basis_set->ndensitymax);

  auto h_rank = Kokkos::create_mirror_view(d_rank);
  auto h_num_ms_combs = Kokkos::create_mirror_view(d_num_ms_combs);
  auto h_offsets = Kokkos::create_mirror_view(d_offsets);
  auto h_mus = Kokkos::create_mirror_view(d_mus);
  auto h_ns = Kokkos::create_mirror_view(d_ns);
  auto h_ls = Kokkos::create_mirror_view(d_ls);
  auto h_ms_combs = Kokkos::create_mirror_view(d_ms_combs);
  auto h_ctildes = Kokkos::create_mirror_view(d_ctildes);

  // copy values on host

  for (int n = 0; n < nelements; n++) {
    const int total_basis_size_rank1 = basis_set->total_basis_size_rank1[n];
    const int total_basis_size = basis_set->total_basis_size[n];

    ACECTildeBasisFunction *basis_rank1 = basis_set->basis_rank1[n];
    ACECTildeBasisFunction *basis = basis_set->basis[n];

    const int ndensity = basis_set->map_embedding_specifications.at(n).ndensity;

    int idx_rho = 0;

    // rank=1
    for (int offset = 0; offset < total_basis_size_rank1; ++offset) {
      ACECTildeBasisFunction *func = &basis_rank1[offset];
      h_rank(n, offset) = 1;
      h_mus(n, offset, 0) = func->mus[0];
      h_ns(n, offset, 0) = func->ns[0];
      for (int p = 0; p < ndensity; p++)
        h_ctildes(n, idx_rho, p) = func->ctildes[p];
      h_offsets(n, idx_rho) = offset;
      idx_rho++;
    }

    // rank > 1
    for (int func_ind = 0; func_ind < total_basis_size; ++func_ind) {
      ACECTildeBasisFunction *func = &basis[func_ind];
      // TODO: check if func->ctildes are zero, then skip

      const int offset = total_basis_size_rank1 + func_ind;

      const int rank = h_rank(n, offset) = func->rank;
      h_num_ms_combs(n, offset) = func->num_ms_combs;
      for (int t = 0; t < rank; t++) {
        h_mus(n, offset, t) = func->mus[t];
        h_ns(n, offset, t) = func->ns[t];
        h_ls(n, offset, t) = func->ls[t];
      }

      // loop over {ms} combinations in sum
      for (int ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind) {
        auto ms = &func->ms_combs[ms_ind * rank]; // current ms-combination (of length = rank)
        for (int t = 0; t < rank; t++)
          h_ms_combs(n, idx_rho, t) = ms[t];

        for (int p = 0; p < ndensity; ++p) {
          // real-part only multiplication
          h_ctildes(n, idx_rho, p) = func->ctildes[ms_ind * ndensity + p];
        }
        h_offsets(n, idx_rho) = offset;
        idx_rho++;
      }
    }
  }

  Kokkos::deep_copy(d_rank, h_rank);
  Kokkos::deep_copy(d_num_ms_combs, h_num_ms_combs);
  Kokkos::deep_copy(d_offsets, h_offsets);
  Kokkos::deep_copy(d_mus, h_mus);
  Kokkos::deep_copy(d_ns, h_ns);
  Kokkos::deep_copy(d_ls, h_ls);
  Kokkos::deep_copy(d_ms_combs, h_ms_combs);
  Kokkos::deep_copy(d_ctildes, h_ctildes);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::init_style()
{
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR,"Pair style pace/kk can currently only run on a single "
                         "CPU thread");

    PairPACE::init_style();
    return;
  }

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style PACE requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style PACE requires newton pair on");

  // neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair pace/kk");

  auto basis_set = aceimpl->basis_set;

  nelements = basis_set->nelements;
  lmax = basis_set->lmax;
  nradmax = basis_set->nradmax;
  nradbase = basis_set->nradbase;

  // spherical harmonics

  MemKK::realloc_kokkos(alm, "pace:alm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(blm, "pace:blm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(cl, "pace:cl", lmax + 1);
  MemKK::realloc_kokkos(dl, "pace:dl", lmax + 1);

  pre_compute_harmonics(lmax);
  copy_pertype();
  copy_splines();
  copy_tilde();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPACEKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairPACE::init_one(i,j);

  k_scale.h_view(i,j) = k_scale.h_view(j,i) = scale[i][j];
  k_scale.template modify<LMPHostType>();

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairPACE::coeff(narg,arg);

  // Set up element lists

  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int i = 1; i <= atom->ntypes; i++)
    h_map(i) = map[i];

  Kokkos::deep_copy(d_map,h_map);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::allocate()
{
  PairPACE::allocate();

  int n = atom->ntypes + 1;
  MemKK::realloc_kokkos(d_map, "pace:map", n);

  MemKK::realloc_kokkos(k_cutsq, "pace:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();

  MemKK::realloc_kokkos(k_scale, "pace:scale", n, n);
  d_scale = k_scale.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& maxneigh) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (maxneigh < num_neighs) maxneigh = num_neighs;
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (host_flag) {
    atomKK->sync(Host,X_MASK|TYPE_MASK);
    PairPACE::compute(eflag_in,vflag_in);
    atomKK->modified(Host,F_MASK);
    return;
  }

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  if (!force->newton_pair)
    error->all(FLERR,"PairPACEKokkos requires 'newton on'");

  if (recursive)
    error->all(FLERR,"Must use 'product' algorithm with pair pace/kk on the GPU");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_scale.template sync<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  maxneigh = 0;
  Kokkos::parallel_reduce("pace::find_maxneigh", inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  int vector_length_default = 1;
  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 32;

  chunk_size = MIN(chunksize,inum); // "chunksize" variable is set by user
  chunk_offset = 0;

  grow(chunk_size, maxneigh);

  EV_FLOAT ev;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    Kokkos::deep_copy(weights, 0.0);
    Kokkos::deep_copy(weights_rank1, 0.0);
    Kokkos::deep_copy(A, 0.0);
    Kokkos::deep_copy(A_rank1, 0.0);
    Kokkos::deep_copy(rhos, 0.0);
    Kokkos::deep_copy(rho_core, 0.0);

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairPACEComputeNeigh>(chunk_size,team_size,vector_length);
      int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh> policy_neigh(chunk_size,team_size,vector_length);
      policy_neigh = policy_neigh.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
    }

    //ComputeRadial
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairPACEComputeRadial>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeRadial> policy_radial(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeRadial",policy_radial,*this);
    }

    //ComputeYlm
    {
      int vector_length = vector_length_default;
      int team_size = 16;
      check_team_size_for<TagPairPACEComputeYlm>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeYlm> policy_ylm(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeYlm",policy_ylm,*this);
    }

    //ComputeAi
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairPACEComputeAi>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi> policy_ai(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeAi",policy_ai,*this);
    }

    //ConjugateAi
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEConjugateAi> policy_conj_ai(0,chunk_size);
      Kokkos::parallel_for("ConjugateAi",policy_conj_ai,*this);
    }

    //ComputeRho
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeRho> policy_rho(0,chunk_size*idx_rho_max);
      Kokkos::parallel_for("ComputeRho",policy_rho,*this);
    }

    //ComputeFS
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeFS> policy_fs(0,chunk_size);
      Kokkos::parallel_for("ComputeFS",policy_fs,*this);
    }

    //ComputeWeights
    {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeWeights> policy_weights(0,chunk_size*idx_rho_max);
      Kokkos::parallel_for("ComputeWeights",policy_weights,*this);
    }

    //ComputeDerivative
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagPairPACEComputeDerivative>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeDerivative> policy_derivative(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeDerivative",policy_derivative,*this);
    }

    //ComputeForce
    {
      if (evflag) {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALF,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALFTHREAD,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce("ComputeForce",policy_force, *this, ev_tmp);
        }
      } else {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALF,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for(policy_force, *this);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALFTHREAD,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for("ComputeForce",policy_force, *this);
        }
      }
    }
    ev += ev_tmp;
    chunk_offset += chunk_size;

  } // end while

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  // get a pointer to scratch memory
  // This is used to cache whether or not an atom is within the cutoff
  // If it is, inside is assigned to 1, otherwise -1
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * maxneigh; // offset into pointer for entire team
  int* inside = (int*)team.team_shmem().get_shmem(team.team_size() * maxneigh * sizeof(int), 0) + scratch_shift;

  // loop over list of all neighbors within force cutoff
  // distsq[] = distance sq to each
  // rlist[] = distance vector to each
  // nearest[] = atom indices of neighbors

  int ncount = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& count) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const int jtype = type(j);

    const F_FLOAT delx = xtmp - x(j,0);
    const F_FLOAT dely = ytmp - x(j,1);
    const F_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    inside[jj] = -1;
    if (rsq < d_cutsq(itype,jtype)) {
     inside[jj] = 1;
     count++;
    }
  },ncount);

  d_ncount(ii) = ncount;

  Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& offset, bool final) {

    if (inside[jj] < 0) return;

    if (final) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const F_FLOAT delx = xtmp - x(j,0);
      const F_FLOAT dely = ytmp - x(j,1);
      const F_FLOAT delz = ztmp - x(j,2);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      const F_FLOAT r = sqrt(rsq);
      const F_FLOAT rinv = 1.0/r;
      const int mu_j = d_map(type(j));
      d_mu(ii,offset) = mu_j;
      d_rnorms(ii,offset) = r;
      d_rhats(ii,offset,0) = -delx*rinv;
      d_rhats(ii,offset,1) = -dely*rinv;
      d_rhats(ii,offset,2) = -delz*rinv;
      d_nearest(ii,offset) = j;
    }
    offset++;
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRadial, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeRadial>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const double r_norm = d_rnorms(ii, jj);
  const int mu_i = d_map(type(i));
  const int mu_j = d_mu(ii, jj);

  evaluate_splines(ii, jj, r_norm, nradbase, nradmax, mu_i, mu_j);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeYlm, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeYlm>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const double xn = d_rhats(ii, jj, 0);
  const double yn = d_rhats(ii, jj, 1);
  const double zn = d_rhats(ii, jj, 2);
  compute_ylm(ii,jj,xn,yn,zn,lmax);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeAi, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const int mu_j = d_mu(ii, jj);

  // rank = 1
  for (int n = 0; n < nradbase; n++)
    Kokkos::atomic_add(&A_rank1(ii, mu_j, n), gr(ii, jj, n) * Y00);

  // rank > 1
  for (int n = 0; n < nradmax; n++) {
    for (int l = 0; l <= lmax; l++) {
      for (int m = 0; m <= l; m++) {
        const int idx = l * (l + 1) + m; // (l, m)
        Kokkos::atomic_add(&A(ii, mu_j, n, idx).re, fr(ii, jj, n, l) * ylm(ii, jj, idx).re);
        Kokkos::atomic_add(&A(ii, mu_j, n, idx).im, fr(ii, jj, n, l) * ylm(ii, jj, idx).im);
      }
    }
  }

  // hard-core repulsion
  Kokkos::atomic_add(&rho_core(ii), cr(ii, jj));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEConjugateAi, const int& ii) const
{
  //complex conjugate A's (for NEGATIVE (-m) terms)
  // for rank > 1
  for (int mu_j = 0; mu_j < nelements; mu_j++) {
    for (int n = 0; n < nradmax; n++) {
      for (int l = 0; l <= lmax; l++) {
        //fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
        for (int m = 1; m <= l; m++) {
          const int idx = l * (l + 1) + m; // (l, m)
          const int idxm = l * (l + 1) - m; // (l, -m)
          const int factor = m % 2 == 0 ? 1 : -1;
          A(ii, mu_j, n, idxm) = A(ii, mu_j, n, idx).conj() * (double)factor;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRho, const int& iter) const
{
  const int idx_rho = iter / chunk_size;
  const int ii = iter % chunk_size;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  if (idx_rho >= d_idx_rho_count(mu_i)) return;

  const int ndensity = d_ndensity(mu_i);

  const int offset = d_offsets(mu_i, idx_rho);
  const int rank = d_rank(mu_i, offset);
  const int r = rank - 1;

  // Basis functions B with iterative product and density rho(p) calculation
  if (rank == 1) {
    const int mu = d_mus(mu_i, offset, 0);
    const int n = d_ns(mu_i, offset, 0);
    double A_cur = A_rank1(ii, mu, n - 1);
    for (int p = 0; p < ndensity; ++p) {
      //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
      Kokkos::atomic_add(&rhos(ii, p), d_ctildes(mu_i, idx_rho, p) * A_cur);
    }
  } else { // rank > 1
    // loop over {ms} combinations in sum

    // loop over m, collect B  = product of A with given ms
    A_forward_prod(ii, idx_rho, 0) = complex::one();

    // fill forward A-product triangle
    for (int t = 0; t < rank; t++) {
      //TODO: optimize ns[t]-1 -> ns[t] during functions construction
      const int mu = d_mus(mu_i, offset, t);
      const int n = d_ns(mu_i, offset, t);
      const int l = d_ls(mu_i, offset, t);
      const int m = d_ms_combs(mu_i, idx_rho, t); // current ms-combination (of length = rank)
      const int idx = l * (l + 1) + m; // (l, m)
      A_list(ii, idx_rho, t) = A(ii, mu, n - 1, idx);
      A_forward_prod(ii, idx_rho, t + 1) = A_forward_prod(ii, idx_rho, t) * A_list(ii, idx_rho, t);
    }

    complex A_backward_prod = complex::one();

    // fill backward A-product triangle
    for (int t = r; t >= 1; t--) {
      const complex dB = A_forward_prod(ii, idx_rho, t) * A_backward_prod; // dB - product of all A's except t-th
      dB_flatten(ii, idx_rho, t) = dB;

      A_backward_prod = A_backward_prod * A_list(ii, idx_rho, t);
    }
    dB_flatten(ii, idx_rho, 0) = A_forward_prod(ii, idx_rho, 0) * A_backward_prod;

    const complex B = A_forward_prod(ii, idx_rho, rank);

    for (int p = 0; p < ndensity; ++p) {
      // real-part only multiplication
      Kokkos::atomic_add(&rhos(ii, p), B.real_part_product(d_ctildes(mu_i, idx_rho, p)));
    }
  }
}

/* ---------------------------------------------------------------------- */


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeFS, const int& ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  const double rho_cut = d_rho_core_cutoff(mu_i);
  const double drho_cut = d_drho_core_cutoff(mu_i);
  const int ndensity = d_ndensity(mu_i);

  double evdwl, fcut, dfcut;
  evdwl = fcut = dfcut = 0.0;

  inner_cutoff(rho_core(ii), rho_cut, drho_cut, fcut, dfcut);
  FS_values_and_derivatives(ii, evdwl, mu_i);

  dF_drho_core(ii) = evdwl * dfcut + 1;
  for (int p = 0; p < ndensity; ++p)
    dF_drho(ii, p) *= fcut;


  // tally energy contribution
  if (eflag) {
    double evdwl_cut = evdwl * fcut + rho_core(ii);
    // E0 shift
    evdwl_cut += d_E0vals(mu_i);
    e_atom(ii) = evdwl_cut;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeWeights, const int& iter) const
{
  const int idx_rho = iter / chunk_size;
  const int ii = iter % chunk_size;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  if (idx_rho >= d_idx_rho_count(mu_i)) return;

  const int ndensity = d_ndensity(mu_i);

  const int offset = d_offsets(mu_i, idx_rho);
  const int rank = d_rank(mu_i, offset);

  // Weights and theta calculation

  if (rank == 1) {
    const int mu = d_mus(mu_i, offset, 0);
    const int n = d_ns(mu_i, offset, 0);
    double theta = 0.0;
    for (int p = 0; p < ndensity; ++p) {
      // for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
      theta += dF_drho(ii, p) * d_ctildes(mu_i, idx_rho, p);
    }
    Kokkos::atomic_add(&weights_rank1(ii, mu, n - 1), theta);
  } else { // rank > 1
    double theta = 0.0;
    for (int p = 0; p < ndensity; ++p)
      theta += dF_drho(ii, p) * d_ctildes(mu_i, idx_rho, p);

    theta *= 0.5; // 0.5 factor due to possible double counting ???
    for (int t = 0; t < rank; ++t) {
      const int m_t = d_ms_combs(mu_i, idx_rho, t);
      const int factor = (m_t % 2 == 0 ? 1 : -1);
      const complex dB = dB_flatten(ii, idx_rho, t);
      const int mu_t = d_mus(mu_i, offset, t);
      const int n_t = d_ns(mu_i, offset, t);
      const int l_t = d_ls(mu_i, offset, t);
      const int idx = l_t * (l_t + 1) + m_t; // (l, m)
      const complex value = theta * dB;
      Kokkos::atomic_add(&(weights(ii, mu_t, n_t - 1, idx).re), value.re);
      Kokkos::atomic_add(&(weights(ii, mu_t, n_t - 1, idx).im), value.im);
      // update -m_t (that could also be positive), because the basis is half_basis
      const int idxm = l_t * (l_t + 1) - m_t; // (l, -m)
      const complex valuem = theta * dB.conj() * (double)factor;
      Kokkos::atomic_add(&(weights(ii, mu_t, n_t - 1, idxm).re), valuem.re);
      Kokkos::atomic_add(&(weights(ii, mu_t, n_t - 1, idxm).im), valuem.im);
    }
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeDerivative, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeDerivative>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int i = d_ilist[ii + chunk_offset];

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  const int itype = type(i);
  const double scale = d_scale(itype,itype);

  const int mu_j = d_mu(ii, jj);
  double r_hat[3];
  r_hat[0] = d_rhats(ii, jj, 0);
  r_hat[1] = d_rhats(ii, jj, 1);
  r_hat[2] = d_rhats(ii, jj, 2);
  const double r = d_rnorms(ii, jj);
  const double rinv = 1.0/r;

  double f_ji[3];
  f_ji[0] = f_ji[1] = f_ji[2] = 0;

  // for rank = 1
  for (int n = 0; n < nradbase; ++n) {
    if (weights_rank1(ii, mu_j, n) == 0) continue;
    double &DG = dgr(ii, jj, n);
    double DGR = DG * Y00;
    DGR *= weights_rank1(ii, mu_j, n);
    f_ji[0] += DGR * r_hat[0];
    f_ji[1] += DGR * r_hat[1];
    f_ji[2] += DGR * r_hat[2];
  }

  // for rank > 1
  for (int n = 0; n < nradmax; n++) {
    for (int l = 0; l <= lmax; l++) {
      const double R_over_r = fr(ii, jj, n, l) * rinv;
      const double DR = dfr(ii, jj, n, l);

      // for m >= 0
      for (int m = 0; m <= l; m++) {
        const int idx = l * (l + 1) + m; // (l, m)
        complex w = weights(ii, mu_j, n, idx);
        if (w.re == 0.0 && w.im == 0.0) continue;
        // counting for -m cases if m > 0
        if (m > 0) {
          w.re *= 2.0;
          w.im *= 2.0;
        }

        complex DY[3];
        DY[0] = dylm(ii, jj, idx, 0);
        DY[1] = dylm(ii, jj, idx, 1);
        DY[2] = dylm(ii, jj, idx, 2);
        const complex Y_DR = ylm(ii, jj, idx) * DR;

        complex grad_phi_nlm[3];
        grad_phi_nlm[0] = Y_DR * r_hat[0] + DY[0] * R_over_r;
        grad_phi_nlm[1] = Y_DR * r_hat[1] + DY[1] * R_over_r;
        grad_phi_nlm[2] = Y_DR * r_hat[2] + DY[2] * R_over_r;
        // real-part multiplication only
        f_ji[0] += w.real_part_product(grad_phi_nlm[0]);
        f_ji[1] += w.real_part_product(grad_phi_nlm[1]);
        f_ji[2] += w.real_part_product(grad_phi_nlm[2]);
      }
    }
  }

  // hard-core repulsion
  const double fpair = dF_drho_core(ii) * dcr(ii,jj);
  f_ij(ii, jj, 0) = scale * f_ji[0] + fpair * r_hat[0];
  f_ij(ii, jj, 1) = scale * f_ji[1] + fpair * r_hat[1];
  f_ij(ii, jj, 2) = scale * f_ji[2] + fpair * r_hat[2];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>, const int& ii, EV_FLOAT& ev) const
{
  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  const auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  const auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type(i);
  const double scale = d_scale(itype,itype);

  const int ncount = d_ncount(ii);

  F_FLOAT fitmp[3] = {0.0,0.0,0.0};
  for (int jj = 0; jj < ncount; jj++) {
    int j = d_nearest(ii,jj);

    double r_hat[3];
    r_hat[0] = d_rhats(ii, jj, 0);
    r_hat[1] = d_rhats(ii, jj, 1);
    r_hat[2] = d_rhats(ii, jj, 2);
    const double r = d_rnorms(ii, jj);
    const double delx = -r_hat[0]*r;
    const double dely = -r_hat[1]*r;
    const double delz = -r_hat[2]*r;

    const double fpairx = f_ij(ii, jj, 0);
    const double fpairy = f_ij(ii, jj, 1);
    const double fpairz = f_ij(ii, jj, 2);

    fitmp[0] += fpairx;
    fitmp[1] += fpairy;
    fitmp[2] += fpairz;
    a_f(j,0) -= fpairx;
    a_f(j,1) -= fpairy;
    a_f(j,2) -= fpairz;

    // tally per-atom virial contribution
    if (EVFLAG && vflag_either)
      v_tally_xyz<NEIGHFLAG>(ev, i, j, fpairx, fpairy, fpairz, delx, dely, delz);
  }

  a_f(i,0) += fitmp[0];
  a_f(i,1) += fitmp[1];
  a_f(i,2) += fitmp[2];

  // tally energy contribution
  if (EVFLAG && eflag_either) {
    const double evdwl = scale*e_atom(ii);
    //ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
    if (eflag_global) ev.evdwl += evdwl;
    if (eflag_atom) d_eatom[i] += evdwl;
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fx, const F_FLOAT &fy, const F_FLOAT &fz,
      const F_FLOAT &delx, const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const E_FLOAT v0 = delx*fx;
  const E_FLOAT v1 = dely*fy;
  const E_FLOAT v2 = delz*fz;
  const E_FLOAT v3 = delx*fy;
  const E_FLOAT v4 = delx*fz;
  const E_FLOAT v5 = dely*fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i,0) += 0.5*v0;
    a_vatom(i,1) += 0.5*v1;
    a_vatom(i,2) += 0.5*v2;
    a_vatom(i,3) += 0.5*v3;
    a_vatom(i,4) += 0.5*v4;
    a_vatom(i,5) += 0.5*v5;
    a_vatom(j,0) += 0.5*v0;
    a_vatom(j,1) += 0.5*v1;
    a_vatom(j,2) += 0.5*v2;
    a_vatom(j,3) += 0.5*v3;
    a_vatom(j,4) += 0.5*v4;
    a_vatom(j,5) += 0.5*v5;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::pre_compute_harmonics(int lmax)
{
  auto h_alm = Kokkos::create_mirror_view(alm);
  auto h_blm = Kokkos::create_mirror_view(blm);
  auto h_cl = Kokkos::create_mirror_view(cl);
  auto h_dl = Kokkos::create_mirror_view(dl);

  for (int l = 1; l <= lmax; l++) {
    const double lsq = l * l;
    const double ld = 2 * l;
    const double l1 = (4 * lsq - 1);
    const double l2 = lsq - ld + 1;
    for (int m = 0; m < l - 1; m++) {
      const double msq = m * m;
      const double a = sqrt((double(l1)) / (double(lsq - msq)));
      const double b = -sqrt((double(l2 - msq)) / (double(4 * l2 - 1)));
      const int idx = l * (l + 1) + m; // (l, m)
      h_alm(idx) = a;
      h_blm(idx) = b;
    }
  }

  for (int l = 1; l <= lmax; l++) {
    h_cl(l) = -sqrt(1.0 + 0.5 / (double(l)));
    h_dl(l) = sqrt(double(2 * (l - 1) + 3));
  }

  Kokkos::deep_copy(alm, h_alm);
  Kokkos::deep_copy(blm, h_blm);
  Kokkos::deep_copy(cl, h_cl);
  Kokkos::deep_copy(dl, h_dl);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_barplm(int ii, int jj, double rz, int lmax) const
{
  // requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!
  // prefactors include 1/sqrt(2) factor compared to reference

  // l=0, m=0
  // plm(ii, jj, 0, 0) = Y00/sq1o4pi; //= sq1o4pi;
  plm(ii, jj, 0) = Y00; //= 1;
  dplm(ii, jj, 0) = 0.0;

  if (lmax > 0) {

    // l=1, m=0
    plm(ii, jj, 2) = Y00 * sq3 * rz;
    dplm(ii, jj, 2) = Y00 * sq3;

    // l=1, m=1
    plm(ii, jj, 3) = -sq3o2 * Y00;
    dplm(ii, jj, 3) = 0.0;

    // loop l = 2, lmax
    for (int l = 2; l <= lmax; l++) {
      for (int m = 0; m < l - 1; m++) {
        const int idx = l * (l + 1) + m; // (l, m)
        const int idx1 = (l - 1) * l + m; // (l - 1, m)
        const int idx2 = (l - 2) * (l - 1) + m; // (l - 2, m)
        plm(ii, jj, idx) = alm(idx) * (rz * plm(ii, jj, idx1) + blm(idx) * plm(ii, jj, idx2));
        dplm(ii, jj, idx) = alm(idx) * (plm(ii, jj, idx1) + rz * dplm(ii, jj, idx1) + blm(idx) * dplm(ii, jj, idx2));
      }
      const int idx = l * (l + 1) + l; // (l, l)
      const int idx1 = l * (l + 1) + l - 1; // (l, l - 1)
      const int idx2 = (l - 1) * l + l - 1; // (l - 1, l - 1)
      const double t = dl(l) * plm(ii, jj, idx2);
      plm(ii, jj, idx1) = t * rz;
      dplm(ii, jj, idx1) = t;
      plm(ii, jj, idx) = cl(l) * plm(ii, jj, idx2);
      dplm(ii, jj, idx) = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_ylm(int ii, int jj, double rx, double ry, double rz, int lmax) const
{
  // requires rx^2 + ry^2 + rz^2 = 1 , NO CHECKING IS PERFORMED !!!!!!!!!

  complex phase;
  complex phasem, mphasem1;
  complex dyx, dyy, dyz;
  complex rdy;

  phase.re = rx;
  phase.im = ry;

  // compute barplm
  compute_barplm(ii, jj, rz, lmax);

  // m = 0
  for (int l = 0; l <= lmax; l++) {
    const int idx = l * (l + 1);

    ylm(ii, jj, idx).re = plm(ii, jj, idx);
    ylm(ii, jj, idx).im = 0.0;

    dyz.re = dplm(ii, jj, idx);
    rdy.re = dyz.re * rz;

    dylm(ii, jj, idx, 0).re = -rdy.re * rx;
    dylm(ii, jj, idx, 0).im = 0.0;
    dylm(ii, jj, idx, 1).re = -rdy.re * ry;
    dylm(ii, jj, idx, 1).im = 0.0;
    dylm(ii, jj, idx, 2).re = dyz.re - rdy.re * rz;
    dylm(ii, jj, idx, 2).im = 0;
  }
  // m = 1
  for (int l = 1; l <= lmax; l++) {
    const int idx = l * (l + 1) + 1;

    ylm(ii, jj, idx) = phase * plm(ii, jj, idx);

    dyx.re = plm(ii, jj, idx);
    dyx.im = 0.0;
    dyy.re = 0.0;
    dyy.im = plm(ii, jj, idx);
    dyz.re = phase.re * dplm(ii, jj, idx);
    dyz.im = phase.im * dplm(ii, jj, idx);

    rdy.re = rx * dyx.re + +rz * dyz.re;
    rdy.im = ry * dyy.im + rz * dyz.im;

    dylm(ii, jj, idx, 0).re = dyx.re - rdy.re * rx;
    dylm(ii, jj, idx, 0).im = -rdy.im * rx;
    dylm(ii, jj, idx, 1).re = -rdy.re * ry;
    dylm(ii, jj, idx, 1).im = dyy.im - rdy.im * ry;
    dylm(ii, jj, idx, 2).re = dyz.re - rdy.re * rz;
    dylm(ii, jj, idx, 2).im = dyz.im - rdy.im * rz;
  }

  // m > 1
  phasem = phase;
  for (int m = 2; m <= lmax; m++) {

    mphasem1.re = phasem.re * double(m);
    mphasem1.im = phasem.im * double(m);
    phasem = phasem * phase;

    for (int l = m; l <= lmax; l++) {
      const int idx = l * (l + 1) + m;

      ylm(ii, jj, idx).re = phasem.re * plm(ii, jj, idx);
      ylm(ii, jj, idx).im = phasem.im * plm(ii, jj, idx);

      dyx = mphasem1 * plm(ii, jj, idx);
      dyy.re = -dyx.im;
      dyy.im = dyx.re;
      dyz = phasem * dplm(ii, jj, idx);

      rdy.re = rx * dyx.re + ry * dyy.re + rz * dyz.re;
      rdy.im = rx * dyx.im + ry * dyy.im + rz * dyz.im;

      dylm(ii, jj, idx, 0).re = dyx.re - rdy.re * rx;
      dylm(ii, jj, idx, 0).im = dyx.im - rdy.im * rx;
      dylm(ii, jj, idx, 1).re = dyy.re - rdy.re * ry;
      dylm(ii, jj, idx, 1).im = dyy.im - rdy.im * ry;
      dylm(ii, jj, idx, 2).re = dyz.re - rdy.re * rz;
      dylm(ii, jj, idx, 2).im = dyz.im - rdy.im * rz;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::cutoff_func_poly(const double r, const double r_in, const double delta_in, double &fc, double &dfc) const
{
  if (r <= r_in-delta_in) {
    fc = 1;
    dfc = 0;
  } else if (r >= r_in ) {
    fc = 0;
    dfc = 0;
  } else {
    double x = 1 - 2 * (1 + (r - r_in) / delta_in);
    fc = 0.5 + 7.5 / 2. * (x / 4. - pow(x, 3) / 6. + pow(x, 5) / 20.);
    dfc = -7.5 / delta_in * (0.25 - x * x / 2.0 + pow(x, 4) / 4.);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::Fexp(const double x, const double m, double &F, double &DF) const
{
  const double w = 1.e6;
  const double eps = 1e-10;

  const double lambda = pow(1.0 / w, m - 1.0);
  if (abs(x) > eps) {
    double g;
    const double a = abs(x);
    const double am = pow(a, m);
    const double w3x3 = pow(w * a, 3); //// use cube
    const double sign_factor = (signbit(x) ? -1 : 1);
    if (w3x3 > 30.0)
        g = 0.0;
    else
        g = exp(-w3x3);

    const double omg = 1.0 - g;
    F = sign_factor * (omg * am + lambda * g * a);
    const double dg = -3.0 * w * w * w * a * a * g;
    DF = m * pow(a, m - 1.0) * omg - am * dg + lambda * dg * a + lambda * g;
  } else {
    F = lambda * x;
    DF = lambda;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::FexpShiftedScaled(const double rho, const double mexp, double &F, double &DF) const
{
  const double eps = 1e-10;

  if (abs(mexp - 1.0) < eps) {
    F = rho;
    DF = 1;
  } else {
    const double a = abs(rho);
    const double exprho = exp(-a);
    const double nx = 1. / mexp;
    const double xoff = pow(nx, (nx / (1.0 - nx))) * exprho;
    const double yoff = pow(nx, (1 / (1.0 - nx))) * exprho;
    const double sign_factor = (signbit(rho) ? -1 : 1);
    F = sign_factor * (pow(xoff + a, mexp) - yoff);
    DF = yoff + mexp * (-xoff + 1.0) * pow(xoff + a, mexp - 1.);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::inner_cutoff(const double rho_core, const double rho_cut, const double drho_cut,
                                     double &fcut, double &dfcut) const
{
  double rho_low = rho_cut - drho_cut;
  if (rho_core >= rho_cut) {
    fcut = 0;
    dfcut = 0;
  } else if (rho_core <= rho_low) {
    fcut = 1;
    dfcut = 0;
  } else {
    cutoff_func_poly(rho_core, rho_cut, drho_cut, fcut, dfcut);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::FS_values_and_derivatives(const int ii, double &evdwl, const int mu_i) const
{
  double F, DF = 0;
  int npoti = d_npoti(mu_i);
  int ndensity = d_ndensity(mu_i);
  for (int p = 0; p < ndensity; p++) {
    const double wpre = d_wpre(mu_i, p);
    const double mexp = d_mexp(mu_i, p);

    if (npoti == FS)
      Fexp(rhos(ii, p), mexp, F, DF);
    else if (npoti == FS_SHIFTEDSCALED)
      FexpShiftedScaled(rhos(ii, p), mexp, F, DF);

    evdwl += F * wpre; // * weight (wpre)
    dF_drho(ii, p) = DF * wpre; // * weight (wpre)
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::evaluate_splines(const int ii, const int jj, double r,
                                                  int /*nradbase_c*/, int /*nradial_c*/,
                                                  int mu_i, int mu_j) const
{
  auto &spline_gk = k_splines_gk.template view<DeviceType>()(mu_i, mu_j);
  auto &spline_rnl = k_splines_rnl.template view<DeviceType>()(mu_i, mu_j);
  auto &spline_hc = k_splines_hc.template view<DeviceType>()(mu_i, mu_j);

  spline_gk.calcSplines(ii, jj, r, gr, dgr);

  spline_rnl.calcSplines(ii, jj, r, d_values, d_derivatives);
  for (int kk = 0; kk < (int)fr.extent(2); kk++) {
    for (int ll = 0; ll < (int)fr.extent(3); ll++) {
      const int flatten = kk*fr.extent(3) + ll;
      fr(ii, jj, kk, ll) = d_values(ii, jj, flatten);
      dfr(ii, jj, kk, ll) = d_derivatives(ii, jj, flatten);
    }
  }

  spline_hc.calcSplines(ii, jj, r, d_values, d_derivatives);
  cr(ii, jj) = d_values(ii, jj, 0);
  dcr(ii, jj) = d_derivatives(ii, jj, 0);
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairPACEKokkos<DeviceType>::SplineInterpolatorKokkos::operator=(const SplineInterpolator &spline) {
    cutoff = spline.cutoff;
    deltaSplineBins = spline.deltaSplineBins;
    ntot = spline.ntot;
    nlut = spline.nlut;
    invrscalelookup = spline.invrscalelookup;
    rscalelookup = spline.rscalelookup;
    num_of_functions = spline.num_of_functions;

    lookupTable = t_ace_3d4("lookupTable", ntot+1, num_of_functions);
    auto h_lookupTable = Kokkos::create_mirror_view(lookupTable);
    for (int i = 0; i < ntot+1; i++)
        for (int j = 0; j < num_of_functions; j++)
            for (int k = 0; k < 4; k++)
                h_lookupTable(i, j, k) = spline.lookupTable(i, j, k);
    Kokkos::deep_copy(lookupTable, h_lookupTable);
}
/* ---------------------------------------------------------------------- */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::SplineInterpolatorKokkos::calcSplines(const int ii, const int jj, const double r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const
{
  double wl, wl2, wl3, w2l1, w3l2;
  double c[4];
  double x = r * rscalelookup;
  int nl = static_cast<int>(floor(x));

  if (nl <= 0)
    Kokkos::abort("Encountered very small distance. Stopping.");

  if (nl < nlut) {
    wl = x - double(nl);
    wl2 = wl * wl;
    wl3 = wl2 * wl;
    w2l1 = 2.0 * wl;
    w3l2 = 3.0 * wl2;
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      for (int idx = 0; idx < 4; idx++)
        c[idx] = lookupTable(nl, func_id, idx);
      d_values(ii, jj, func_id) = c[0] + c[1] * wl + c[2] * wl2 + c[3] * wl3;
      d_derivatives(ii, jj, func_id) = (c[1] + c[2] * w2l1 + c[3] * w3l2) * rscalelookup;
    }
  } else { // fill with zeroes
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      d_values(ii, jj, func_id) = 0.0;
      d_derivatives(ii, jj, func_id) = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairPACEKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairPACEKokkos<DeviceType>::check_team_size_reduce(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType>
template<typename scratch_type>
int PairPACEKokkos<DeviceType>::scratch_size_helper(int values_per_team) {
  typedef Kokkos::View<scratch_type*, Kokkos::DefaultExecutionSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPACEKokkos<DeviceType>::memory_usage()
{
  double bytes = 0;

  bytes += MemKK::memory_usage(A);
  bytes += MemKK::memory_usage(A_rank1);
  bytes += MemKK::memory_usage(A_list);
  bytes += MemKK::memory_usage(A_forward_prod);
  bytes += MemKK::memory_usage(e_atom);
  bytes += MemKK::memory_usage(rhos);
  bytes += MemKK::memory_usage(dF_drho);
  bytes += MemKK::memory_usage(weights);
  bytes += MemKK::memory_usage(weights_rank1);
  bytes += MemKK::memory_usage(rho_core);
  bytes += MemKK::memory_usage(dF_drho_core);
  bytes += MemKK::memory_usage(dB_flatten);
  bytes += MemKK::memory_usage(fr);
  bytes += MemKK::memory_usage(dfr);
  bytes += MemKK::memory_usage(gr);
  bytes += MemKK::memory_usage(dgr);
  bytes += MemKK::memory_usage(d_values);
  bytes += MemKK::memory_usage(d_derivatives);
  bytes += MemKK::memory_usage(cr);
  bytes += MemKK::memory_usage(dcr);
  bytes += MemKK::memory_usage(plm);
  bytes += MemKK::memory_usage(dplm);
  bytes += MemKK::memory_usage(ylm);
  bytes += MemKK::memory_usage(dylm);
  bytes += MemKK::memory_usage(d_ncount);
  bytes += MemKK::memory_usage(d_mu);
  bytes += MemKK::memory_usage(d_rhats);
  bytes += MemKK::memory_usage(d_rnorms);
  bytes += MemKK::memory_usage(d_nearest);
  bytes += MemKK::memory_usage(f_ij);
  bytes += MemKK::memory_usage(d_rho_core_cutoff);
  bytes += MemKK::memory_usage(d_drho_core_cutoff);
  bytes += MemKK::memory_usage(d_E0vals);
  bytes += MemKK::memory_usage(d_ndensity);
  bytes += MemKK::memory_usage(d_npoti);
  bytes += MemKK::memory_usage(d_wpre);
  bytes += MemKK::memory_usage(d_mexp);
  bytes += MemKK::memory_usage(d_idx_rho_count);
  bytes += MemKK::memory_usage(d_rank);
  bytes += MemKK::memory_usage(d_num_ms_combs);
  bytes += MemKK::memory_usage(d_offsets);
  bytes += MemKK::memory_usage(d_mus);
  bytes += MemKK::memory_usage(d_ns);
  bytes += MemKK::memory_usage(d_ls);
  bytes += MemKK::memory_usage(d_ms_combs);
  bytes += MemKK::memory_usage(d_ctildes);
  bytes += MemKK::memory_usage(alm);
  bytes += MemKK::memory_usage(blm);
  bytes += MemKK::memory_usage(cl);
  bytes += MemKK::memory_usage(dl);

  if (k_splines_gk.h_view.data()) {
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        bytes += k_splines_gk.h_view(i, j).memory_usage();
        bytes += k_splines_rnl.h_view(i, j).memory_usage();
        bytes += k_splines_hc.h_view(i, j).memory_usage();
      }
    }
  }

  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairPACEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairPACEKokkos<LMPHostType>;
#endif
}
