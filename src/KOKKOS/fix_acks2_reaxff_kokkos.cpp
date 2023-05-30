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
------------------------------------------------------------------------- */

#include "fix_acks2_reaxff_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_reaxff_kokkos.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 0.0001
#define EV_TO_KCAL_PER_MOL 14.4

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

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | Q_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = Q_MASK | X_MASK;

  nmax = m_cap = 0;
  allocated_flag = 0;
  nprev = 4;

  memory->destroy(s_hist);
  memory->destroy(s_hist_X);
  memory->destroy(s_hist_last);
  grow_arrays(atom->nmax);
  memoryKK->create_kokkos(k_s_hist_last,s_hist_last,2,nprev,"acks2/reax:s_hist_last");
  d_s_hist_last = k_s_hist_last.template view<DeviceType>();
  buf = new double[2*nprev];
  prev_last_rows_rank = 0;

  d_mfill_offset = typename AT::t_int_scalar("acks2/kk:mfill_offset");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixACKS2ReaxFFKokkos<DeviceType>::~FixACKS2ReaxFFKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_s_hist,s_hist);
  memoryKK->destroy_kokkos(k_s_hist_X,s_hist_X);
  memoryKK->destroy_kokkos(k_s_hist_last,s_hist_last);
  delete [] buf;

  deallocate_array();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init()
{
  atomKK->k_q.modify<LMPHostType>();
  atomKK->k_q.sync<DeviceType>();

  FixACKS2ReaxFF::init();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag_qeq;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();

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
  init_hist();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init_shielding_k()
{
  int i,j;
  int ntypes = atom->ntypes;

  k_shield = DAT::tdual_ffloat_2d("acks2/kk:shield",ntypes+1,ntypes+1);
  d_shield = k_shield.template view<DeviceType>();

  for( i = 1; i <= ntypes; ++i )
    for( j = 1; j <= ntypes; ++j )
      k_shield.h_view(i,j) = pow( gamma[i] * gamma[j], -1.5 );

  k_shield.template modify<LMPHostType>();
  k_shield.template sync<DeviceType>();

  k_bcut = DAT::tdual_ffloat_2d("acks2/kk:bcut",ntypes+1,ntypes+1);
  d_bcut = k_bcut.template view<DeviceType>();

  for( i = 1; i <= ntypes; ++i )
    for( j = 1; j <= ntypes; ++j )
      k_bcut.h_view(i,j) = 0.5*(bcut_acks2[i] + bcut_acks2[j]);

  k_bcut.template modify<LMPHostType>();
  k_bcut.template sync<DeviceType>();

  k_tap = DAT::tdual_ffloat_1d("acks2/kk:tap",8);
  d_tap = k_tap.template view<DeviceType>();

  for (i = 0; i < 8; i ++)
    k_tap.h_view(i) = Tap[i];

  k_tap.template modify<LMPHostType>();
  k_tap.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::init_hist()
{
  k_s_hist.clear_sync_state();
  k_s_hist_X.clear_sync_state();
  k_s_hist_last.clear_sync_state();

  Kokkos::deep_copy(d_s_hist,0.0);
  Kokkos::deep_copy(d_s_hist_X,0.0);
  Kokkos::deep_copy(d_s_hist_last,0.0);

  k_s_hist.template modify<DeviceType>();
  k_s_hist_X.template modify<DeviceType>();
  k_s_hist_last.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  atomKK->sync(execution_space,datamask_read);

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
  k_bcut.template sync<DeviceType>();
  k_tap.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  nn = list->inum;
  NN = atom->nlocal + atom->nghost;

  copymode = 1;

  // allocate

  allocate_array();

  if (!allocated_flag || last_allocate < neighbor->lastcall) {

    // get max number of neighbor

    allocate_matrix();

    // last_rows_rank proc must not own zero atoms (unless no atoms total)
    //  otherwise some loops are no-ops and last rows contribution won't
    //  be computed correctly

    int flag = comm->me;
    if (nn == 0) flag = MAXSMALLINT;
    MPI_Allreduce(&flag, &last_rows_rank, 1, MPI_INT, MPI_MIN, world);
    last_rows_flag = (comm->me == last_rows_rank);

    // pass along "s" array history if necessary

    if (prev_last_rows_rank != last_rows_rank) {

      MPI_Request request;
      if (comm->me == last_rows_rank)
        MPI_Irecv(buf,2*nprev,MPI_DOUBLE,
                  prev_last_rows_rank,0,world,&request);

      if (comm->me == prev_last_rows_rank) {

        // pack buffer
        k_s_hist_last.template sync<LMPHostType>();
        auto h_s_hist_last = k_s_hist_last.h_view;
        int n = 0;
        for (int k = 0; k < nprev; k++) {
          buf[n++] = h_s_hist_last(0,k);
          buf[n++] = h_s_hist_last(1,k);
        }

        MPI_Send(buf,2*nprev,MPI_DOUBLE,last_rows_rank,0,world);
      }

      if (comm->me == last_rows_rank) {
        MPI_Wait(&request,MPI_STATUS_IGNORE);

        // unpack buffer
        k_s_hist_last.template sync<LMPHostType>();
        auto h_s_hist_last = k_s_hist_last.h_view;
        int n = 0;
        for (int k = 0; k < nprev; k++) {
          h_s_hist_last(0,k) = buf[n++];
          h_s_hist_last(1,k) = buf[n++];
        }
        k_s_hist_last.template modify<LMPHostType>();
      }
    }

    prev_last_rows_rank = last_rows_rank;

    last_allocate = update->ntimestep;
  }

  // compute_H

  if (execution_space == Host) { // CPU
    if (neighflag == FULL) {
      FixACKS2ReaxFFKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(this);
      Kokkos::parallel_scan(nn,computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixACKS2ReaxFFKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(this);
      Kokkos::parallel_scan(nn,computeH_functor);
    }
  } else { // GPU, use teams
    Kokkos::deep_copy(d_mfill_offset,0);

    int vector_length = 32;
    int atoms_per_team = 4;
    int num_teams = nn / atoms_per_team + (nn % atoms_per_team ? 1 : 0);

    Kokkos::TeamPolicy<DeviceType> policy(num_teams, atoms_per_team,
                                          vector_length);

    if (neighflag == FULL) {
      FixACKS2ReaxFFKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixACKS2ReaxFFKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    }
  }

  need_dup = lmp->kokkos->need_dup<DeviceType>(1);

  if (need_dup)
    dup_X_diag = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated> (d_X_diag); // allocate duplicated memory
  else
    ndup_X_diag = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated> (d_X_diag);

  // compute_X

  Kokkos::deep_copy(d_X_diag,0.0);

  if (execution_space == Host || 1) { // CPU
    if (neighflag == FULL) {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, FULL> computeX_functor(this);
      Kokkos::parallel_scan(nn,computeX_functor);
    } else if (neighflag == HALFTHREAD) {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, HALFTHREAD> computeX_functor(this);
      Kokkos::parallel_scan(nn,computeX_functor);
    } else {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, HALF> computeX_functor(this);
      Kokkos::parallel_scan(nn,computeX_functor);
    }
  } else { // GPU, use teams
    Kokkos::deep_copy(d_mfill_offset,0);

    int vector_length = 32;
    int atoms_per_team = 4;
    int num_teams = nn / atoms_per_team + (nn % atoms_per_team ? 1 : 0);

    Kokkos::TeamPolicy<DeviceType> policy(num_teams, atoms_per_team,
                                          vector_length);
    if (neighflag == FULL) {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, FULL> computeX_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeX_functor);
    } else if (neighflag == HALFTHREAD) {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, HALFTHREAD> computeX_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeX_functor);
    } else {
      FixACKS2ReaxFFKokkosComputeXFunctor<DeviceType, HALF> computeX_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeX_functor);
    }
  }

  if (need_dup) {
    Kokkos::Experimental::contribute(d_X_diag, dup_X_diag);

    // free duplicated memory

    dup_X_diag = decltype(dup_X_diag)();
  }

  if (neighflag != FULL) {
    pack_flag = 4;
    //comm->reverse_comm(this); //Coll_Vector( X_diag );
    k_X_diag.template modify<DeviceType>();
    k_X_diag.template sync<LMPHostType>();
    comm->reverse_comm(this);
    k_X_diag.template modify<LMPHostType>();
    k_X_diag.template sync<DeviceType>();
  }

  if (efield) get_chi_field();

  // init_matvec

  k_s_hist.template sync<DeviceType>();
  k_s_hist_X.template sync<DeviceType>();
  k_s_hist_last.template sync<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2InitMatvec>(0,nn),*this);

  pack_flag = 2;
  // comm->forward_comm(this); //Dist_vector( s );
  k_s.template modify<DeviceType>();
  k_s.template sync<LMPHostType>();
  comm->forward_comm(this);
  more_forward_comm(k_s.h_view.data());
  k_s.template modify<LMPHostType>();
  k_s.template sync<DeviceType>();

  // bicgstab solve over b_s, s

  bicgstab_solve();

  calculate_Q();

  k_s_hist.template modify<DeviceType>();
  k_s_hist_X.template modify<DeviceType>();
  k_s_hist_last.template modify<DeviceType>();

  copymode = 0;

  if (!allocated_flag)
    allocated_flag = 1;

  atomKK->modified(execution_space,datamask_modify);
  k_s.modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::num_neigh_item(int ii, int &maxneigh) const
{
  const int i = d_ilist[ii];
  maxneigh += d_numneigh[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::allocate_matrix()
{
  nmax = atom->nmax;

  // determine the total space for the H matrix

  m_cap = 0;

  // limit scope of functor to allow deallocation of views
  {
    FixACKS2ReaxFFKokkosNumNeighFunctor<DeviceType> neigh_functor(this);
    Kokkos::parallel_reduce(nn,neigh_functor,m_cap);
  }

  // deallocate first to reduce memory overhead

  d_firstnbr = typename AT::t_int_1d();
  d_numnbrs = typename AT::t_int_1d();
  d_jlist = typename AT::t_int_1d();
  d_val = typename AT::t_ffloat_1d();

  d_firstnbr_X = typename AT::t_int_1d();
  d_numnbrs_X = typename AT::t_int_1d();
  d_jlist_X = typename AT::t_int_1d();
  d_val_X = typename AT::t_ffloat_1d();

  // H matrix

  d_firstnbr = typename AT::t_int_1d("acks2/kk:firstnbr",nmax);
  d_numnbrs = typename AT::t_int_1d("acks2/kk:numnbrs",nmax);
  d_jlist = typename AT::t_int_1d("acks2/kk:jlist",m_cap);
  d_val = typename AT::t_ffloat_1d("acks2/kk:val",m_cap);

  // X matrix

  d_firstnbr_X = typename AT::t_int_1d("acks2/kk:firstnbr_X",nmax);
  d_numnbrs_X = typename AT::t_int_1d("acks2/kk:numnbrs_X",nmax);
  d_jlist_X = typename AT::t_int_1d("acks2/kk:jlist_X",m_cap);
  d_val_X = typename AT::t_ffloat_1d("acks2/kk:val_X",m_cap);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::allocate_array()
{
  // 0 to nn-1: owned atoms related to H matrix
  // nn to NN-1: ghost atoms related to H matrix
  // NN to NN+nn-1: owned atoms related to X matrix
  // NN+nn to 2*NN-1: ghost atoms related X matrix
  // 2*NN to 2*NN+1: last two rows, owned by proc 0

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    int size = nmax*2 + 2;

    d_q = typename AT::t_ffloat_1d("acks2/kk:q",size);

    memoryKK->create_kokkos(k_s,s,size,"acks2/kk:s");
    d_s = k_s.template view<DeviceType>();

    d_b_s = typename AT::t_ffloat_1d("acks2/kk:b_s",size);

    d_Hdia_inv = typename AT::t_ffloat_1d("acks2/kk:Hdia_inv",nmax);

    memoryKK->create_kokkos(k_chi_field,chi_field,nmax,"acks2/kk:chi_field");
    d_chi_field = k_chi_field.template view<DeviceType>();

    memoryKK->create_kokkos(k_X_diag,X_diag,nmax,"acks2/kk:X_diag");
    d_X_diag = k_X_diag.template view<DeviceType>();

    d_Xdia_inv = typename AT::t_ffloat_1d("acks2/kk:Xdia_inv",nmax);

    d_p = typename AT::t_ffloat_1d("acks2/kk:p",size);
    d_r = typename AT::t_ffloat_1d("acks2/kk:r",size);

    memoryKK->create_kokkos(k_d,d,size,"acks2/kk:d");
    d_d = k_d.template view<DeviceType>();

    d_g = typename AT::t_ffloat_1d("acks2/kk:g",size);

    memoryKK->create_kokkos(k_q_hat,q_hat,size,"acks2/kk:q_hat");
    d_q_hat = k_q_hat.template view<DeviceType>();

    d_r_hat = typename AT::t_ffloat_1d("acks2/kk:r_hat",size);

    memoryKK->create_kokkos(k_y,y,size,"acks2/kk:y");
    d_y = k_y.template view<DeviceType>();

    memoryKK->create_kokkos(k_z,z,size,"acks2/kk:z");
    d_z = k_z.template view<DeviceType>();
  }

  if (efield) get_chi_field();

  // init_storage
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2Zero>(0,nn),*this);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::deallocate_array()
{
  memoryKK->destroy_kokkos(k_s,s);
  memoryKK->destroy_kokkos(k_chi_field,chi_field);
  memoryKK->destroy_kokkos(X_diag);
  memoryKK->destroy_kokkos(k_d,d);
  memoryKK->destroy_kokkos(k_q_hat,q_hat);
  memoryKK->destroy_kokkos(k_y,y);
  memoryKK->destroy_kokkos(k_z,z);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Zero, const int &ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_s[i] = -params(itype).chi - d_chi_field[i];
    d_s[i] = 0.0;
    d_p[i] = 0.0;
    d_r[i] = 0.0;
    d_d[i] = 0.0;
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::compute_h_item(int ii, int &m_fill, const bool &final) const
{
  const int i = d_ilist[ii];
  int j,jj,jtype;

  if (mask[i] & groupbit) {

    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int itype = type(i);
    const tagint itag = tag(i);
    const int jnum = d_numneigh[i];
    if (final)
      d_firstnbr[i] = m_fill;

    for (jj = 0; jj < jnum; jj++) {
      j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      jtype = type(j);

      const X_FLOAT delx = x(j,0) - xtmp;
      const X_FLOAT dely = x(j,1) - ytmp;
      const X_FLOAT delz = x(j,2) - ztmp;

      if (NEIGHFLAG != FULL) {
        // skip half of the interactions
        const tagint jtag = tag(j);
        if (j >= nlocal) {
          if (itag > jtag) {
            if ((itag+jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag+jtag) % 2 == 1) continue;
          } else {
            if (x(j,2) < ztmp) continue;
            if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
            if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
          }
        }
      }

      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq > cutsq) continue;

      if (final) {
        const F_FLOAT r = sqrt(rsq);
        d_jlist(m_fill) = j;
        const F_FLOAT shldij = d_shield(itype,jtype);
        d_val(m_fill) = calculate_H_k(r,shldij);
      }
      m_fill++;
    }
    if (final)
      d_numnbrs[i] = m_fill - d_firstnbr[i];
  }
}

/* ---------------------------------------------------------------------- */

// Calculate Qeq matrix H where H is a sparse matrix and H[i][j] represents the electrostatic interaction coefficients on atom-i with atom-j
// d_val     - contains the non-zero entries of sparse matrix H
// d_numnbrs - d_numnbrs[i] contains the # of non-zero entries in the i-th row of H (which also represents the # of neighbor atoms with electrostatic interaction coefficients with atom-i)
// d_firstnbr- d_firstnbr[i] contains the beginning index from where the H matrix entries corresponding to row-i is stored in d_val
// d_jlist   - contains the column index corresponding to each entry in d_val

template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::compute_h_team(
    const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
    int atoms_per_team, int vector_length) const {

  // scratch space setup
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_ilist(team.team_shmem(), atoms_per_team);
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_numnbrs(team.team_shmem(), atoms_per_team);
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_firstnbr(team.team_shmem(), atoms_per_team);

  Kokkos::View<int **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_jtype(team.team_shmem(), atoms_per_team, vector_length);
  Kokkos::View<int **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_jlist(team.team_shmem(), atoms_per_team, vector_length);
  Kokkos::View<F_FLOAT **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_r(team.team_shmem(), atoms_per_team, vector_length);

  // team of threads work on atoms with index in [firstatom, lastatom)
  int firstatom = team.league_rank() * atoms_per_team;
  int lastatom =
      (firstatom + atoms_per_team < nn) ? (firstatom + atoms_per_team) : nn;

  // kokkos-thread-0 is used to load info from global memory into scratch space
  if (team.team_rank() == 0) {

    // copy atom indices from d_ilist[firstatom:lastatom] to scratch space s_ilist[0:atoms_per_team]
    // copy # of neighbor atoms for all the atoms with indices in d_ilist[firstatom:lastatom] from d_numneigh to scratch space s_numneigh[0:atoms_per_team]
    // calculate total number of neighbor atoms for all atoms assigned to the current team of threads (Note - Total # of neighbor atoms here provides the
    // upper bound space requirement to store the H matrix values corresponding to the atoms with indices in d_ilist[firstatom:lastatom])

    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team, atoms_per_team),
                          [&](const int &idx, int &totalnbrs, bool final) {
                            int ii = firstatom + idx;

                            if (ii < nn) {
                              const int i = d_ilist[ii];
                              int jnum = d_numneigh[i];

                              if (final) {
                                s_ilist[idx] = i;
                                s_numnbrs[idx] = jnum;
                                s_firstnbr[idx] = totalnbrs;
                              }
                              totalnbrs += jnum;
                            } else {
                              s_numnbrs[idx] = 0;
                            }
                          });
  }

  // barrier ensures that the data moved to scratch space is visible to all the
  // threads of the corresponding team
  team.team_barrier();

  // calculate the global memory offset from where the H matrix values to be
  // calculated by the current team will be stored in d_val
  int team_firstnbr_idx = 0;
  Kokkos::single(Kokkos::PerTeam(team),
                 [=](int &val) {
                   int totalnbrs = s_firstnbr[lastatom - firstatom - 1] +
                                   s_numnbrs[lastatom - firstatom - 1];
                   val = Kokkos::atomic_fetch_add(&d_mfill_offset(), totalnbrs);
                 },
                 team_firstnbr_idx);

  // map the H matrix computation of each atom to kokkos-thread (one atom per
  // kokkos-thread) neighbor computation for each atom is assigned to vector
  // lanes of the corresponding thread
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, atoms_per_team), [&](const int &idx) {
        int ii = firstatom + idx;

        if (ii < nn) {
          const int i = s_ilist[idx];

          if (mask[i] & groupbit) {
            const X_FLOAT xtmp = x(i, 0);
            const X_FLOAT ytmp = x(i, 1);
            const X_FLOAT ztmp = x(i, 2);
            const int itype = type(i);
            tagint itag = tag(i);
            int jnum = s_numnbrs[idx];

            // calculate the write-offset for atom-i's first neighbor
            int atomi_firstnbr_idx = team_firstnbr_idx + s_firstnbr[idx];
            Kokkos::single(Kokkos::PerThread(team),
                           [&]() { d_firstnbr[i] = atomi_firstnbr_idx; });

            // current # of neighbor atoms with non-zero electrostatic
            // interaction coefficients with atom-i which represents the # of
            // non-zero elements in row-i of H matrix
            int atomi_nbrs_inH = 0;

            // calculate H matrix values corresponding to atom-i where neighbors
            // are processed in batches and the batch size is vector_length
            for (int jj_start = 0; jj_start < jnum; jj_start += vector_length) {

              int atomi_nbr_writeIdx = atomi_firstnbr_idx + atomi_nbrs_inH;

              // count the # of neighbor atoms with non-zero electrostatic
              // interaction coefficients with atom-i in the current batch
              int atomi_nbrs_curbatch = 0;

              // compute rsq, jtype, j and store in scratch space which is
              // reused later
              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(team, vector_length),
                  [&](const int &idx, int &m_fill) {
                    const int jj = jj_start + idx;

                    // initialize: -1 represents no interaction with atom-j
                    // where j = d_neighbors(i,jj)
                    s_jlist(team.team_rank(), idx) = -1;

                    if (jj < jnum) {
                      int j = d_neighbors(i, jj);
                      j &= NEIGHMASK;
                      const int jtype = type(j);

                      const X_FLOAT delx = x(j, 0) - xtmp;
                      const X_FLOAT dely = x(j, 1) - ytmp;
                      const X_FLOAT delz = x(j, 2) - ztmp;

                      // valid nbr interaction
                      bool valid = true;
                      if (NEIGHFLAG != FULL) {
                        // skip half of the interactions
                        const tagint jtag = tag(j);
                        if (j >= nlocal) {
                          if (itag > jtag) {
                            if ((itag + jtag) % 2 == 0)
                              valid = false;
                          } else if (itag < jtag) {
                            if ((itag + jtag) % 2 == 1)
                              valid = false;
                          } else {
                            if (x(j, 2) < ztmp)
                              valid = false;
                            if (x(j, 2) == ztmp && x(j, 1) < ytmp)
                              valid = false;
                            if (x(j, 2) == ztmp && x(j, 1) == ytmp && x(j, 0) < xtmp)
                              valid = false;
                          }
                        }
                      }

                      const F_FLOAT rsq =
                          delx * delx + dely * dely + delz * delz;
                      if (rsq > cutsq)
                        valid = false;

                      if (valid) {
                        s_jlist(team.team_rank(), idx) = j;
                        s_jtype(team.team_rank(), idx) = jtype;
                        s_r(team.team_rank(), idx) = sqrt(rsq);
                        m_fill++;
                      }
                    }
                  },
                  atomi_nbrs_curbatch);

              // write non-zero entries of H to global memory
              Kokkos::parallel_scan(
                  Kokkos::ThreadVectorRange(team, vector_length),
                  [&](const int &idx, int &m_fill, bool final) {
                    int j = s_jlist(team.team_rank(), idx);
                    if (final) {
                      if (j != -1) {
                        const int jtype = s_jtype(team.team_rank(), idx);
                        const F_FLOAT r = s_r(team.team_rank(), idx);
                        const F_FLOAT shldij = d_shield(itype, jtype);

                        d_jlist[atomi_nbr_writeIdx + m_fill] = j;
                        d_val[atomi_nbr_writeIdx + m_fill] =
                            calculate_H_k(r, shldij);
                      }
                    }

                    if (j != -1) {
                      m_fill++;
                    }
                  });
              atomi_nbrs_inH += atomi_nbrs_curbatch;
            }

            Kokkos::single(Kokkos::PerThread(team),
                           [&]() { d_numnbrs[i] = atomi_nbrs_inH; });
          }
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
  denom = pow(denom,1.0/3.0);

  return taper * EV_TO_KCAL_PER_MOL / denom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::compute_x_item(int ii, int &m_fill, const bool &final) const
{
  // The X_diag array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_X_diag = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_X_diag),decltype(ndup_X_diag)>::get(dup_X_diag,ndup_X_diag);
  auto a_X_diag = v_X_diag.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  int j,jj,jtype;
  F_FLOAT tmp = 0.0;

  if (mask[i] & groupbit) {

    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int itype = type(i);
    const tagint itag = tag(i);
    const int jnum = d_numneigh[i];
    if (final)
      d_firstnbr_X[i] = m_fill;

    for (jj = 0; jj < jnum; jj++) {
      j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      jtype = type(j);

      const X_FLOAT delx = x(j,0) - xtmp;
      const X_FLOAT dely = x(j,1) - ytmp;
      const X_FLOAT delz = x(j,2) - ztmp;

      if (NEIGHFLAG != FULL) {
        // skip half of the interactions
        const tagint jtag = tag(j);
        if (j >= nlocal) {
          if (itag > jtag) {
            if ((itag+jtag) % 2 == 0) continue;
          } else if (itag < jtag) {
            if ((itag+jtag) % 2 == 1) continue;
          } else {
            if (x(j,2) < ztmp) continue;
            if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
            if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
          }
        }
      }

      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq > cutsq) continue;

      const F_FLOAT bcutoff = d_bcut(itype,jtype);
      const F_FLOAT bcutoff2 = bcutoff*bcutoff;
      if (rsq > bcutoff2) continue;

      if (final) {
        const F_FLOAT r = sqrt(rsq);
        d_jlist_X(m_fill) = j;
        const F_FLOAT X_val = calculate_X_k(r,bcutoff);
        d_val_X(m_fill) = X_val;
        tmp -= X_val;
        if (NEIGHFLAG != FULL)
          a_X_diag[j] -= X_val;
      }
      m_fill++;
    }
    if (final) {
      a_X_diag[i] += tmp;
      d_numnbrs_X[i] = m_fill - d_firstnbr_X[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::compute_x_team(
    const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
    int atoms_per_team, int vector_length) const {

  // The X_diag array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_X_diag = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_X_diag),decltype(ndup_X_diag)>::get(dup_X_diag,ndup_X_diag);
  auto a_X_diag = v_X_diag.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // scratch space setup
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_ilist(team.team_shmem(), atoms_per_team);
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_numnbrs(team.team_shmem(), atoms_per_team);
  Kokkos::View<int *, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_firstnbr(team.team_shmem(), atoms_per_team);

  Kokkos::View<int **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_jtype(team.team_shmem(), atoms_per_team, vector_length);
  Kokkos::View<int **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_jlist(team.team_shmem(), atoms_per_team, vector_length);
  Kokkos::View<F_FLOAT **, Kokkos::ScratchMemorySpace<DeviceType>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      s_r(team.team_shmem(), atoms_per_team, vector_length);

  // team of threads work on atoms with index in [firstatom, lastatom)
  int firstatom = team.league_rank() * atoms_per_team;
  int lastatom =
      (firstatom + atoms_per_team < nn) ? (firstatom + atoms_per_team) : nn;

  // kokkos-thread-0 is used to load info from global memory into scratch space
  if (team.team_rank() == 0) {

   // copy atom indices from d_ilist[firstatom:lastatom] to scratch space s_ilist[0:atoms_per_team]
    // copy # of neighbor atoms for all the atoms with indices in d_ilist[firstatom:lastatom] from d_numneigh to scratch space s_numneigh[0:atoms_per_team]
    // calculate total number of neighbor atoms for all atoms assigned to the current team of threads (Note - Total # of neighbor atoms here provides the
    // upper bound space requirement to store the H matrix values corresponding to the atoms with indices in d_ilist[firstatom:lastatom])

    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team, atoms_per_team),
                          [&](const int &idx, int &totalnbrs, bool final) {
                            int ii = firstatom + idx;

                            if (ii < nn) {
                              const int i = d_ilist[ii];
                              int jnum = d_numneigh[i];

                              if (final) {
                                s_ilist[idx] = i;
                                s_numnbrs[idx] = jnum;
                                s_firstnbr[idx] = totalnbrs;
                              }
                              totalnbrs += jnum;
                            } else {
                              s_numnbrs[idx] = 0;
                            }
                          });
  }

  // barrier ensures that the data moved to scratch space is visible to all the
  // threads of the corresponding team
  team.team_barrier();

  // calculate the global memory offset from where the H matrix values to be
  // calculated by the current team will be stored in d_val_X
  int team_firstnbr_idx = 0;
  Kokkos::single(Kokkos::PerTeam(team),
                 [=](int &val) {
                   int totalnbrs = s_firstnbr[lastatom - firstatom - 1] +
                                   s_numnbrs[lastatom - firstatom - 1];
                   val = Kokkos::atomic_fetch_add(&d_mfill_offset(), totalnbrs);
                 },
                 team_firstnbr_idx);

  // map the H matrix computation of each atom to kokkos-thread (one atom per
  // kokkos-thread) neighbor computation for each atom is assigned to vector
  // lanes of the corresponding thread
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, atoms_per_team), [&](const int &idx) {
        int ii = firstatom + idx;

        if (ii < nn) {
          const int i = s_ilist[idx];

          if (mask[i] & groupbit) {
            const X_FLOAT xtmp = x(i, 0);
            const X_FLOAT ytmp = x(i, 1);
            const X_FLOAT ztmp = x(i, 2);
            const int itype = type(i);
            tagint itag = tag(i);
            int jnum = s_numnbrs[idx];

            // calculate the write-offset for atom-i's first neighbor
            int atomi_firstnbr_idx = team_firstnbr_idx + s_firstnbr[idx];
            Kokkos::single(Kokkos::PerThread(team),
                           [&]() { d_firstnbr_X[i] = atomi_firstnbr_idx; });

            // current # of neighbor atoms with non-zero electrostatic
            // interaction coefficients with atom-i which represents the # of
            // non-zero elements in row-i of H matrix
            int atomi_nbrs_inH = 0;

            // calculate H matrix values corresponding to atom-i where neighbors
            // are processed in batches and the batch size is vector_length
            for (int jj_start = 0; jj_start < jnum; jj_start += vector_length) {

              int atomi_nbr_writeIdx = atomi_firstnbr_idx + atomi_nbrs_inH;

              // count the # of neighbor atoms with non-zero electrostatic
              // interaction coefficients with atom-i in the current batch
              int atomi_nbrs_curbatch = 0;

              // compute rsq, jtype, j and store in scratch space which is
              // reused later
              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(team, vector_length),
                  [&](const int &idx, int &m_fill) {
                    const int jj = jj_start + idx;

                    // initialize: -1 represents no interaction with atom-j
                    // where j = d_neighbors(i,jj)
                    s_jlist(team.team_rank(), idx) = -1;

                    if (jj < jnum) {
                      int j = d_neighbors(i, jj);
                      j &= NEIGHMASK;
                      const int jtype = type(j);

                      const X_FLOAT delx = x(j, 0) - xtmp;
                      const X_FLOAT dely = x(j, 1) - ytmp;
                      const X_FLOAT delz = x(j, 2) - ztmp;

                      // valid nbr interaction
                      bool valid = true;
                      if (NEIGHFLAG != FULL) {
                        // skip half of the interactions
                        const tagint jtag = tag(j);
                        if (j >= nlocal) {
                          if (itag > jtag) {
                            if ((itag + jtag) % 2 == 0)
                              valid = false;
                          } else if (itag < jtag) {
                            if ((itag + jtag) % 2 == 1)
                              valid = false;
                          } else {
                            if (x(j, 2) < ztmp)
                              valid = false;
                            if (x(j, 2) == ztmp && x(j, 1) < ytmp)
                              valid = false;
                            if (x(j, 2) == ztmp && x(j, 1) == ytmp &&
                              x(j, 0) < xtmp)
                              valid = false;
                          }
                        }
                      }

                      const F_FLOAT rsq =
                          delx * delx + dely * dely + delz * delz;
                      if (rsq > cutsq)
                        valid = false;

                       const F_FLOAT bcutoff = d_bcut(itype,jtype);
                       const F_FLOAT bcutoff2 = bcutoff*bcutoff;
                       if (rsq > bcutoff2)
                         valid = false;

                      if (valid) {
                        s_jlist(team.team_rank(), idx) = j;
                        s_jtype(team.team_rank(), idx) = jtype;
                        s_r(team.team_rank(), idx) = sqrt(rsq);
                        m_fill++;
                      }
                    }
                  },
                  atomi_nbrs_curbatch);

              // write non-zero entries of H to global memory
              Kokkos::parallel_scan(
                  Kokkos::ThreadVectorRange(team, vector_length),
                  [&](const int &idx, int &m_fill, bool final) {
                    int j = s_jlist(team.team_rank(), idx);
                    if (final) {
                      if (j != -1) {
                        const int jtype = s_jtype(team.team_rank(), idx);
                        const F_FLOAT r = s_r(team.team_rank(), idx);
                        const F_FLOAT bcutoff = d_bcut(itype, jtype);

                        d_jlist_X[atomi_nbr_writeIdx + m_fill] = j;
                        const F_FLOAT X_val = calculate_X_k(r, bcutoff);
                        d_val_X[atomi_nbr_writeIdx + m_fill] =
                            X_val;
                        a_X_diag[i] -= X_val;
                        if (NEIGHFLAG != FULL)
                          a_X_diag[j] -= X_val;
                      }
                    }

                    if (j != -1) {
                      m_fill++;
                    }
                  });
              atomi_nbrs_inH += atomi_nbrs_curbatch;
            }

            Kokkos::single(Kokkos::PerThread(team),
                           [&]() { d_numnbrs_X[i] = atomi_nbrs_inH; });
          }
        }
      });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixACKS2ReaxFFKokkos<DeviceType>::calculate_X_k( const double &r, const double &bcut) const
{
  const F_FLOAT d = r/bcut;
  const F_FLOAT d3 = d*d*d;
  const F_FLOAT omd = 1.0 - d;
  const F_FLOAT omd2 = omd*omd;
  const F_FLOAT omd6 = omd2*omd2*omd2;

  return bond_softness*d3*omd6;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2InitMatvec, const int &ii) const
{
  if (d_X_diag[ii] == 0.0)
    d_Xdia_inv[ii] = 1.0;
  else
    d_Xdia_inv[ii] = 1.0 / d_X_diag[ii];

  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_s[i] = -params(itype).chi - d_chi_field[i];
    d_b_s[NN+i] = 0.0;

    d_s[i] = 4*(d_s_hist(i,0)+d_s_hist(i,2))-(6*d_s_hist(i,1)+d_s_hist(i,3));
    d_s[NN+i] = 4*(d_s_hist_X(i,0)+d_s_hist_X(i,2))-(6*d_s_hist_X(i,1)+d_s_hist_X(i,3));
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    for (int k = 0; k < 2; k++) {
      d_b_s[2*NN+k] = 0.0;
      d_s[2*NN+k] = 4*(d_s_hist_last(k,0)+d_s_hist_last(k,2))-(6*d_s_hist_last(k,1)+d_s_hist_last(k,3));
    }
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixACKS2ReaxFFKokkos<DeviceType>::bicgstab_solve()
{
  int i;
  F_FLOAT my_norm,norm_sqr,my_dot,dot_sqr;
  double tmp, sigma, rho, rho_old, rnorm, bnorm;

  // sparse_matvec( &H, &X, x, d );
  sparse_matvec_acks2(d_s, d_d);

  pack_flag = 1;
  k_d.template modify<DeviceType>();
  k_d.template sync<LMPHostType>();
  if (neighflag != FULL)
    comm->reverse_comm(this); //Coll_vector( d );
  more_reverse_comm(k_d.h_view.data());
  k_d.template modify<LMPHostType>();
  k_d.template sync<DeviceType>();

  // vector_sum( r , 1.,  b, -1., d, nn );
  // bnorm = parallel_norm( b, nn );
  my_norm = 0.0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Norm1>(0,nn),*this,my_norm);
  norm_sqr = 0.0;
  MPI_Allreduce( &my_norm, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  bnorm = sqrt(norm_sqr);

  // rnorm = parallel_norm( r, nn);
  my_norm = 0.0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Norm2>(0,nn),*this,my_norm);
  norm_sqr = 0.0;
  MPI_Allreduce( &my_norm, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  rnorm = sqrt(norm_sqr);

  if (bnorm == 0.0 ) bnorm = 1.0;
  Kokkos::deep_copy(d_r_hat,d_r);
  omega = 1.0;
  rho = 1.0;

  for (i = 1; i < imax && rnorm / bnorm > tolerance; ++i) {
    // rho = parallel_dot( r_hat, r, nn);
    my_dot = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Dot1>(0,nn),*this,my_dot);
    dot_sqr = 0.0;
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    rho = dot_sqr;
    if (rho == 0.0) break;

    if (i > 1) {
      beta = (rho / rho_old) * (alpha / omega);

      // vector_sum( p , 1., r, beta, q, nn);
      // vector_sum( q , 1., p, -omega, z, nn);
      // pre-conditioning
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2Precon1A>(0,nn),*this);
    } else {

      // vector_copy( p , r nn);
      // pre-conditioning
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2Precon1B>(0,nn),*this);
    }

    pack_flag = 1;
    // comm->forward_comm(this); //Dist_vector( d );
    k_d.template modify<DeviceType>();
    k_d.template sync<LMPHostType>();
    comm->forward_comm(this);
    more_forward_comm(k_d.h_view.data());
    k_d.template modify<LMPHostType>();
    k_d.template sync<DeviceType>();

    // sparse_matvec( &H, &X, d, z );
    sparse_matvec_acks2(d_d, d_z);

    pack_flag = 2;
    k_z.template modify<DeviceType>();
    k_z.template sync<LMPHostType>();
    if (neighflag != FULL)
      comm->reverse_comm(this); //Coll_vector( z );
    more_reverse_comm(k_z.h_view.data());
    k_z.template modify<LMPHostType>();
    k_z.template sync<DeviceType>();

    // tmp = parallel_dot( r_hat, z, nn);
    my_dot = dot_sqr = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Dot2>(0,nn),*this,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    tmp = dot_sqr;
    alpha = rho / tmp;

    // vector_sum( q, 1., r, -alpha, z, nn);
    // tmp = parallel_dot( q, q, nn);
    my_dot = dot_sqr = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Dot3>(0,nn),*this,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    tmp = dot_sqr;

    // early convergence check
    if (tmp < tolerance) {
      // vector_add( x, alpha, d, nn);
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2Add>(0,nn),*this);
      break;
    }

    // pre-conditioning
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2Precon2>(0,nn),*this);

    // sparse_matvec( &H, &X, q_hat, y );
    pack_flag = 3;
    // comm->forward_comm(this); //Dist_vector( q_hat );
    k_q_hat.template modify<DeviceType>();
    k_q_hat.template sync<LMPHostType>();
    comm->forward_comm(this);
    more_forward_comm(k_q_hat.h_view.data());
    k_q_hat.template modify<LMPHostType>();
    k_q_hat.template sync<DeviceType>();

    sparse_matvec_acks2(d_q_hat, d_y);

    pack_flag = 3;
    k_y.template modify<DeviceType>();
    k_y.template sync<LMPHostType>();
    if (neighflag != FULL)
      comm->reverse_comm(this); //Coll_vector( y );
    more_reverse_comm(k_y.h_view.data());
    k_y.template modify<LMPHostType>();
    k_y.template sync<DeviceType>();

    // sigma = parallel_dot( y, q, nn);
    my_dot = dot_sqr = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Dot4>(0,nn),*this,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    sigma = dot_sqr;

    // tmp = parallel_dot( y, y, nn);
    my_dot = dot_sqr = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Dot5>(0,nn),*this,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    tmp = dot_sqr;

    omega = sigma / tmp;

    // vector_sum( g , alpha, d, omega, q_hat, nn);
    // vector_add( x, 1., g, nn);
    // vector_sum( r , 1., q, -omega, y, nn);
    // rnorm = parallel_norm( r, nn);
    my_dot = dot_sqr = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagACKS2Norm3>(0,nn),*this,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    rnorm = sqrt(dot_sqr);

    if (omega == 0) break;
    rho_old = rho;
  }

  if (comm->me == 0) {
    if (omega == 0 || rho == 0) {
      error->warning(FLERR,"Fix acks2/reaxff/kk BiCGStab numerical breakdown, omega = {:.8}, rho = {:.8}",
                      omega,rho);
    } else if (i >= imax) {
      error->warning(FLERR,"Fix acks2/reaxff/kk BiCGStab convergence failed after {} iterations "
                           "at step {}", i, update->ntimestep);
    }
  }

  return i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::calculate_Q()
{
  pack_flag = 2;
  //comm->forward_comm( this ); //Dist_vector( s );
  k_s.modify<DeviceType>();
  k_s.sync<LMPHostType>();
  comm->forward_comm(this);
  k_s.modify<LMPHostType>();
  k_s.sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2CalculateQ>(0,NN),*this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::sparse_matvec_acks2(typename AT::t_ffloat_1d &d_xx_in, typename AT::t_ffloat_1d &d_bb_in)
{
  d_xx = d_xx_in;
  d_bb = d_bb_in;

  if (need_dup)
    dup_bb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated> (d_bb); // allocate duplicated memory
  else
    ndup_bb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated> (d_bb);

  Kokkos::deep_copy(d_bb,0.0); // can make more efficient?
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2SparseMatvec1>(0,nn),*this);

  if (neighflag == FULL) {
    int teamsize;
    if (execution_space == Host) teamsize = 1;
    else teamsize = 128;

    Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceType,TagACKS2SparseMatvec3_Full>(nn,teamsize),*this);
  } else if (neighflag == HALFTHREAD)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2SparseMatvec3_Half<HALFTHREAD> >(0,nn),*this);
  else if (neighflag == HALF)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagACKS2SparseMatvec3_Half<HALF> >(0,nn),*this);

  if (need_dup) {
    Kokkos::Experimental::contribute(d_bb, dup_bb);

    // free duplicated memory

    dup_bb = decltype(dup_bb)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2SparseMatvec1, const int &ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);
  if (mask[i] & groupbit) {
    d_bb[i] = params(itype).eta * d_xx[i];
    d_bb[NN + i] = d_X_diag[i] * d_xx[NN + i];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2SparseMatvec2, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_bb[i] = 0.0;
    d_bb[NN + i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2SparseMatvec3_Half<NEIGHFLAG>, const int &ii) const
{
  // The bb array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_bb = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_bb),decltype(ndup_bb)>::get(dup_bb,ndup_bb);
  auto a_bb = v_bb.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    F_FLOAT tmp = 0.0;

    // H Matrix
    for(int jj = d_firstnbr[i]; jj < d_firstnbr[i] + d_numnbrs[i]; jj++) {
      const int j = d_jlist(jj);
      tmp += d_val(jj) * d_xx[j];
      a_bb[j] += d_val(jj) * d_xx[i];
    }
    a_bb[i] += tmp;

    // X Matrix
    tmp = 0.0;
    for(int jj = d_firstnbr_X[i]; jj < d_firstnbr_X[i] + d_numnbrs_X[i]; jj++) {
      const int j = d_jlist_X(jj);
      tmp += d_val_X(jj) * d_xx[NN + j];
     a_bb[NN + j] += d_val_X(jj) * d_xx[NN + i];
    }
    a_bb[NN + i] += tmp;

    // Identity Matrix
    a_bb[NN + i] += d_xx[i];
    a_bb[i] += d_xx[NN + i];

    // Second-to-last row/column
    a_bb[2*NN] += d_xx[NN + i];
    a_bb[NN + i] += d_xx[2*NN];

    // Last row/column
    a_bb[2*NN + 1] += d_xx[i];
    a_bb[i] += d_xx[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2SparseMatvec3_Full, const membertype &team) const
{
  int i = d_ilist[team.league_rank()];
  if (mask[i] & groupbit) {
    F_FLOAT sum;
    F_FLOAT sum2;

    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT &sum) {
      const int j = d_jlist(jj);
      sum += d_val(jj) * d_xx[j];
    }, sum);
    team.team_barrier();

    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, d_firstnbr_X[i], d_firstnbr_X[i] + d_numnbrs_X[i]), [&] (const int &jj, F_FLOAT &sum2) {
      const int j = d_jlist_X(jj);
      sum2 += d_val_X(jj) * d_xx[NN + j];
    }, sum2);

    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      d_bb[i] += sum;
      d_bb[NN + i] += sum2;

      // Identity Matrix
      d_bb[NN + i] += d_xx[i];
      d_bb[i] += d_xx[NN + i];

      // Second-to-last row/column
      Kokkos::atomic_add(&(d_bb[2*NN]),d_xx[NN + i]);
      d_bb[NN + i] += d_xx[2*NN];

      // Last row/column
      Kokkos::atomic_add(&(d_bb[2*NN + 1]),d_xx[i]);
      d_bb[i] += d_xx[2*NN + 1];
    });
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Norm1, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_r[i] = d_b_s[i] - d_d[i];
    d_r[NN+i] = d_b_s[NN+i] - d_d[NN+i];

    lsum += d_b_s[i] * d_b_s[i];
    lsum += d_b_s[NN+i] * d_b_s[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_r[2*NN] = d_b_s[2*NN] - d_d[2*NN];
    d_r[2*NN + 1] = d_b_s[2*NN + 1] - d_d[2*NN + 1];

    lsum += d_b_s[2*NN] * d_b_s[2*NN];
    lsum += d_b_s[2*NN + 1] * d_b_s[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Norm2, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    lsum += d_r[i] * d_r[i];
    lsum += d_r[NN + i] * d_r[NN + i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    lsum += d_r[2*NN] * d_r[2*NN];
    lsum += d_r[2*NN + 1] * d_r[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Dot1, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    lsum += d_r_hat[i] * d_r[i];
    lsum += d_r_hat[NN+i] * d_r[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    lsum += d_r_hat[2*NN] * d_r[2*NN];
    lsum += d_r_hat[2*NN + 1] * d_r[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Precon1A, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
      d_q[i] = d_p[i] - omega*d_z[i];
      d_q[NN+i] = d_p[NN+i] - omega*d_z[NN+i];

      d_p[i] = d_r[i] + beta*d_q[i];
      d_p[NN+i] = d_r[NN+i] + beta*d_q[NN+i];

      d_d[i] = d_p[i]*d_Hdia_inv[i];
      d_d[NN+i] = d_p[NN+i]*d_Xdia_inv[i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_q[2*NN] = d_p[2*NN] - omega*d_z[2*NN];
    d_q[2*NN + 1] = d_p[2*NN + 1] - omega*d_z[2*NN + 1];

    d_p[2*NN] = d_r[2*NN] + beta*d_q[2*NN];
    d_p[2*NN + 1] = d_r[2*NN + 1] + beta*d_q[2*NN + 1];

    d_d[2*NN] = d_p[2*NN];
    d_d[2*NN + 1] = d_p[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Precon1B, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
      d_p[i] = d_r[i] ;
      d_p[NN+i] = d_r[NN+i];

      d_d[i] = d_p[i]*d_Hdia_inv[i];
      d_d[NN+i] = d_p[NN+i]*d_Xdia_inv[i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_p[2*NN] = d_r[2*NN];
    d_p[2*NN + 1] = d_r[2*NN + 1];

    d_d[2*NN] = d_p[2*NN];
    d_d[2*NN + 1] = d_p[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Dot2, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    lsum += d_r_hat[i] * d_z[i];
    lsum += d_r_hat[NN+i] * d_z[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    lsum += d_r_hat[2*NN] * d_z[2*NN];
    lsum += d_r_hat[2*NN + 1] * d_z[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Dot3, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_q[i] = d_r[i] - alpha*d_z[i];
    d_q[NN+i] = d_r[NN+i] - alpha*d_z[NN+i];

    lsum += d_q[i] * d_q[i];
    lsum += d_q[NN+i] * d_q[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_q[2*NN] = d_r[2*NN] - alpha*d_z[2*NN];
    d_q[2*NN + 1] = d_r[2*NN + 1] - alpha*d_z[2*NN + 1];

    lsum += d_q[2*NN] * d_q[2*NN];
    lsum += d_q[2*NN + 1] * d_q[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Dot4, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    lsum += d_y[i] * d_q[i];
    lsum += d_y[NN+i] * d_q[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    lsum += d_y[2*NN] * d_q[2*NN];
    lsum += d_y[2*NN + 1] * d_q[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Dot5, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    lsum += d_y[i] * d_y[i];
    lsum += d_y[NN+i] * d_y[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    lsum += d_y[2*NN] * d_y[2*NN];
    lsum += d_y[2*NN + 1] * d_y[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Add, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_s[i] += alpha * d_d[i];
    d_s[NN+i] += alpha * d_d[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_s[2*NN] += alpha*d_d[2*NN];
    d_s[2*NN + 1] += alpha*d_d[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Precon2, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_q_hat[i] = d_q[i]*d_Hdia_inv[i];
    d_q_hat[NN+i] = d_q[NN+i]*d_Xdia_inv[i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_q_hat[2*NN] = d_q[2*NN];
    d_q_hat[2*NN + 1] = d_q[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2Norm3, const int &ii, double &lsum) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_g[i] = alpha*d_d[i] + omega*d_q_hat[i];
    d_g[NN+i] = alpha*d_d[NN+i] + omega*d_q_hat[NN+i];

    d_s[i] += d_g[i];
    d_s[NN+i] += d_g[NN+i];

    d_r[i] = d_q[i] - omega*d_y[i];
    d_r[NN+i] = d_q[NN+i] - omega*d_y[NN+i];

    lsum += d_r[i] * d_r[i];
    lsum += d_r[NN+i] * d_r[NN+i];
  }

  // last two rows
  if (last_rows_flag && ii == 0) {
    d_g[2*NN] = alpha*d_d[2*NN] + omega*d_q_hat[2*NN];
    d_g[2*NN + 1] = alpha*d_d[2*NN + 1] + omega*d_q_hat[2*NN + 1];

    d_s[2*NN] += d_g[2*NN];
    d_s[2*NN + 1] += d_g[2*NN + 1];

    d_r[2*NN] = d_q[2*NN] - omega*d_y[2*NN];
    d_r[2*NN + 1] = d_q[2*NN + 1] - omega*d_y[2*NN + 1];

    lsum += d_r[2*NN] * d_r[2*NN];
    lsum += d_r[2*NN + 1] * d_r[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixACKS2ReaxFFKokkos<DeviceType>::operator() (TagACKS2CalculateQ, const int &i) const
{
  if (mask[i] & groupbit) {

    q(i) = d_s(i);

    if (i < nlocal) {
      /* backup s */
      for (int k = nprev-1; k > 0; --k) {
        d_s_hist(i,k) = d_s_hist(i,k-1);
        d_s_hist_X(i,k) = d_s_hist_X(i,k-1);
      }
      d_s_hist(i,0) = d_s[i];
      d_s_hist_X(i,0) = d_s[NN+i];
    }
  }

  // last two rows
  if (last_rows_flag && i == 0) {
    for (int i = 0; i < 2; ++i) {
      for (int k = nprev-1; k > 0; --k)
        d_s_hist_last(i,k) = d_s_hist_last(i,k-1);
      d_s_hist_last(i,0) = d_s[2*NN+i];
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixACKS2ReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes;

  int size = 2*nmax + 2;

  bytes = size*nprev * sizeof(double); // s_hist
  bytes += nmax*4 * sizeof(double); // storage
  bytes += size*11 * sizeof(double); // storage
  bytes += n_cap*4 * sizeof(int); // matrix...
  bytes += m_cap*2 * sizeof(int);
  bytes += m_cap*2 * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_s_hist.template sync<LMPHostType>();
  k_s_hist_X.template sync<LMPHostType>();

  k_s_hist.template modify<LMPHostType>(); // force reallocation on host
  k_s_hist_X.template modify<LMPHostType>();

  memoryKK->grow_kokkos(k_s_hist,s_hist,nmax,nprev,"acks2:s_hist");
  memoryKK->grow_kokkos(k_s_hist_X,s_hist_X,nmax,nprev,"acks2:s_hist_X");

  d_s_hist = k_s_hist.template view<DeviceType>();
  d_s_hist_X = k_s_hist_X.template view<DeviceType>();

  k_s_hist.template modify<LMPHostType>();
  k_s_hist_X.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
{
  k_s_hist.template sync<LMPHostType>();
  k_s_hist_X.template sync<LMPHostType>();

  FixACKS2ReaxFF::copy_arrays(i,j,delflag);

  k_s_hist.template modify<LMPHostType>();
  k_s_hist_X.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_s_hist.sync_device();
  k_s_hist_X.sync_device();

  Sorter.sort(LMPDeviceType(), k_s_hist.d_view);
  Sorter.sort(LMPDeviceType(), k_s_hist_X.d_view);

  k_s_hist.modify_device();
  k_s_hist_X.modify_device();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixACKS2ReaxFFKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_s_hist.template sync<LMPHostType>();
  k_s_hist_X.template sync<LMPHostType>();

  return FixACKS2ReaxFF::pack_exchange(i,buf);
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixACKS2ReaxFFKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_s_hist.template sync<LMPHostType>();
  k_s_hist_X.template sync<LMPHostType>();

  int n = FixACKS2ReaxFF::unpack_exchange(nlocal,buf);

  k_s_hist.template modify<LMPHostType>();
  k_s_hist_X.template modify<LMPHostType>();

  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixACKS2ReaxFFKokkos<DeviceType>::get_chi_field()
{
  atomKK->sync(Host,X_MASK|MASK_MASK|IMAGE_MASK);
  FixQEqReaxFF::get_chi_field();
  k_chi_field.modify_host();
  k_chi_field.sync_device();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixACKS2ReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixACKS2ReaxFFKokkos<LMPHostType>;
#endif
}
