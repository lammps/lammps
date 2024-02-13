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
   Contributing authors: Ray Shan (SNL), Stan Moore (SNL),
                          Kamesh Arumugam (NVIDIA)

   Nicholas Curtis (AMD), Leopold Grinberd (AMD), and Gina Sitaraman (AMD):
     - Reduced math overhead: enabled specialized calls (e.g., cbrt for a
         cube root instead of pow) and use power/exponential laws to reduce the
         number of exponentials evaluated, etc.
     - Fused the CG solve for "S" and "T" matrices
     - Improved the SpMV algorithm by using vector instead of team level
         parallelism on GPUs
------------------------------------------------------------------------- */

#include "fix_qeq_reaxff_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec_kokkos.h"
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
FixQEqReaxFFKokkos<DeviceType>::
FixQEqReaxFFKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixQEqReaxFF(lmp, narg, arg)
{
  kokkosable = 1;
  comm_forward = comm_reverse = 2; // fused
  forward_comm_device = exchange_comm_device = sort_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | Q_MASK | MASK_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = X_MASK;

  nmax = m_cap = 0;
  allocated_flag = 0;
  nprev = 4;
  maxexchange = nprev*2;

  memory->destroy(s_hist);
  memory->destroy(t_hist);
  grow_arrays(atom->nmax);

  d_mfill_offset = typename AT::t_int_scalar("qeq/kk:mfill_offset");

  converged = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxFFKokkos<DeviceType>::~FixQEqReaxFFKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_s_hist,s_hist);
  memoryKK->destroy_kokkos(k_t_hist,t_hist);
  memoryKK->destroy_kokkos(k_chi_field,chi_field);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init()
{
  atomKK->sync(execution_space,Q_MASK);

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
    ("FixQEqReaxFF::params",ntypes+1);
  params = k_params.template view<DeviceType>();

  for (int n = 1; n <= ntypes; n++) {
    k_params.h_view(n).chi = chi[n];
    k_params.h_view(n).eta = eta[n];
    k_params.h_view(n).gamma = gamma[n];
  }
  k_params.template modify<LMPHostType>();

  cutsq = swb * swb;

  init_shielding_k();
  init_hist();

  last_allocate = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init_shielding_k()
{
  int i,j;
  int ntypes = atom->ntypes;

  k_shield = DAT::tdual_ffloat_2d("qeq/kk:shield",ntypes+1,ntypes+1);
  d_shield = k_shield.template view<DeviceType>();

  for (i = 1; i <= ntypes; ++i)
    for (j = 1; j <= ntypes; ++j)
      k_shield.h_view(i,j) = pow(gamma[i] * gamma[j], -1.5);

  k_shield.template modify<LMPHostType>();
  k_shield.template sync<DeviceType>();

  k_tap = DAT::tdual_ffloat_1d("qeq/kk:tap",8);
  d_tap = k_tap.template view<DeviceType>();

  for (i = 0; i < 8; i ++)
    k_tap.h_view(i) = Tap[i];

  k_tap.template modify<LMPHostType>();
  k_tap.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::init_hist()
{
  k_s_hist.clear_sync_state();
  k_t_hist.clear_sync_state();

  Kokkos::deep_copy(d_s_hist,0.0);
  Kokkos::deep_copy(d_t_hist,0.0);

  k_s_hist.template modify<DeviceType>();
  k_t_hist.template modify<DeviceType>();
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
  k_tap.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  nn = list->inum;

  copymode = 1;

  // allocate

  allocate_array();

  // get max number of neighbor

  if (!allocated_flag || last_allocate < neighbor->lastcall) {
    allocate_matrix();
    last_allocate = update->ntimestep;
  }

  // compute_H

  if (execution_space == Host) { // CPU
    if (neighflag == FULL) {
      FixQEqReaxFFKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(this);
      Kokkos::parallel_scan(nn,computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixQEqReaxFFKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(this);
      Kokkos::parallel_scan(nn,computeH_functor);
    }
  } else { // GPU, use teams
    Kokkos::deep_copy(d_mfill_offset,0);

    int atoms_per_team = FixQEqReaxFFKokkos<DeviceType>::compute_h_teamsize;
    int vector_length = FixQEqReaxFFKokkos<DeviceType>::compute_h_vectorsize;

    int num_teams = nn / atoms_per_team + (nn % atoms_per_team ? 1 : 0);

    Kokkos::TeamPolicy<DeviceType> policy(num_teams, atoms_per_team,
                                          vector_length);
    if (neighflag == FULL) {
      FixQEqReaxFFKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixQEqReaxFFKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    }
  }

  // init_matvec

  if (efield) get_chi_field();

  k_s_hist.template sync<DeviceType>();
  k_t_hist.template sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqInitMatvec>(0,nn),*this);

  pack_flag = 2;
  k_st.template modify<DeviceType>();
  comm->forward_comm(this);
  k_st.template sync<DeviceType>();

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_o = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(d_o); // allocate duplicated memory
  else
    ndup_o = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(d_o);

  //  cg solve over b_s, s & b_t, t

  matvecs = cg_solve();

  // calculate_Q();

  k_s_hist.template sync<DeviceType>();
  k_t_hist.template sync<DeviceType>();
  calculate_q();
  k_s_hist.template modify<DeviceType>();
  k_t_hist.template modify<DeviceType>();

  copymode = 0;

  if (!allocated_flag)
    allocated_flag = 1;

  // free duplicated memory

  if (need_dup)
    dup_o = decltype(dup_o)();

  atomKK->modified(execution_space,datamask_modify);

  d_neighbors = typename AT::t_neighbors_2d();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::num_neigh_item(int ii, int &maxneigh) const
{
  const int i = d_ilist[ii];
  maxneigh += d_numneigh[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::allocate_matrix()
{
  nmax = atom->nmax;

  // determine the total space for the H matrix

  m_cap = 0;

  // limit scope of functor to allow deallocation of views
  {
    FixQEqReaxFFKokkosNumNeighFunctor<DeviceType> neigh_functor(this);
    Kokkos::parallel_reduce(nn,neigh_functor,m_cap);
  }

  // deallocate first to reduce memory overhead

  d_firstnbr = typename AT::t_int_1d();
  d_numnbrs = typename AT::t_int_1d();
  d_jlist = typename AT::t_int_1d();
  d_val = typename AT::t_ffloat_1d();

  d_firstnbr = typename AT::t_int_1d("qeq/kk:firstnbr",nmax);
  d_numnbrs = typename AT::t_int_1d("qeq/kk:numnbrs",nmax);
  d_jlist = typename AT::t_int_1d("qeq/kk:jlist",m_cap);
  d_val = typename AT::t_ffloat_1d("qeq/kk:val",m_cap);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::allocate_array()
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;

    k_o = DAT::tdual_ffloat2_1d("qeq/kk:o",nmax);
    d_o = k_o.template view<DeviceType>();
    h_o = k_o.h_view;

    d_p = typename AT::t_ffloat2_1d("qeq/kk:p",nmax);
    d_r = typename AT::t_ffloat2_1d("qeq/kk:r",nmax);
    k_d = DAT::tdual_ffloat2_1d("qeq/kk:d",nmax);
    d_d = k_d.template view<DeviceType>();
    h_d = k_d.h_view;

    d_Hdia_inv = typename AT::t_ffloat_1d("qeq/kk:Hdia_inv",nmax);

    d_b_st = typename AT::t_ffloat2_1d("qeq/kk:b_st",nmax);

    k_st = DAT::tdual_ffloat2_1d("qeq/kk:st",nmax);
    d_st = k_st.template view<DeviceType>();
    h_st = k_st.h_view;

    memoryKK->create_kokkos(k_chi_field,chi_field,nmax,"qeq/kk:chi_field");
    d_chi_field = k_chi_field.template view<DeviceType>();
  }

  // init_storage

  if (efield) get_chi_field();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqZero>(0,nn),*this);

}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqZero, const int &ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_st(i,0) = -params(itype).chi - d_chi_field[i];
    d_b_st(i,1) = -1.0;
    d_st(i,0) = 0.0;
    d_st(i,1) = 0.0;
    d_p(i,0) = 0.0;
    d_p(i,1) = 0.0;
    d_o(i,0) = 0.0;
    d_o(i,1) = 0.0;
    d_r(i,0) = 0.0;
    d_r(i,1) = 0.0;
    d_d(i,0) = 0.0;
    d_d(i,1) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::compute_h_item(int ii, int &m_fill, const bool &final) const
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
void FixQEqReaxFFKokkos<DeviceType>::compute_h_team(
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
            tagint itag = tag(i); // removed "const" to work around GCC 7 bug
            int jnum = s_numnbrs[idx]; // removed "const" to work around GCC 7 bug

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
double FixQEqReaxFFKokkos<DeviceType>::calculate_H_k(const F_FLOAT &r, const F_FLOAT &shld) const
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
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqInitMatvec, const int &ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_st(i,0) = -params(itype).chi - d_chi_field[i];
    d_b_st(i,1) = -1.0;
    d_st(i,0) = 4*(d_s_hist(i,0)+d_s_hist(i,2))-(6*d_s_hist(i,1)+d_s_hist(i,3));
    d_st(i,1) = d_t_hist(i,2) + 3*(d_t_hist(i,0) - d_t_hist(i,1));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::cg_solve()
// b = b_st, x = st;
{
  converged = 0;

  F_FLOAT2 tmp;
  F_FLOAT2 sig_old;
  F_FLOAT2 b_norm;

  // sparse_matvec(&H, x, q);
  sparse_matvec_kokkos(d_st);

  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    comm->reverse_comm(this); //Coll_vector(q);
    k_o.template sync<DeviceType>();
  }

  // vector_sum(r , 1.,  b, -1., q, nn);
  // preconditioning: d[j] = r[j] * Hdia_inv[j];
  // b_norm = parallel_norm(b, nn);
  F_FLOAT2 my_norm;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagQEqNorm1>(0,nn),*this,my_norm);
  F_FLOAT2 norm_sqr;
  MPI_Allreduce(&my_norm.v, &norm_sqr.v, 2, MPI_DOUBLE, MPI_SUM, world);
  b_norm.v[0] = sqrt(norm_sqr.v[0]);
  b_norm.v[1] = sqrt(norm_sqr.v[1]);

  // sig_new = parallel_dot(r, d, nn);
  F_FLOAT2 my_dot;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagQEqDot1>(0,nn),*this,my_dot);
  F_FLOAT2 dot_sqr;
  MPI_Allreduce(&my_dot.v, &dot_sqr.v, 2, MPI_DOUBLE, MPI_SUM, world);
  F_FLOAT2 sig_new;
  sig_new = dot_sqr;

  F_FLOAT residual[2] = {0.0, 0.0};
  int loop;
  for (loop = 1; (loop < imax); loop++) {
    if (!(converged & 1))
      residual[0] = sqrt(sig_new.v[0]) / b_norm.v[0];
    if (!(converged & 2))
      residual[1] = sqrt(sig_new.v[1]) / b_norm.v[1];
    converged = static_cast<int>(residual[0] <= tolerance) | (static_cast<int>(residual[1] <= tolerance) << 1);

    if (converged == 3) {
      // both cg solves have converged
      break;
    }

    // comm->forward_comm(this); //Dist_vector(d);
    pack_flag = 1;
    // mark size 2 for
    k_d.template modify<DeviceType>();
    comm->forward_comm(this, 2);
    k_d.template sync<DeviceType>();

    // sparse_matvec(&H, d, q);
    sparse_matvec_kokkos(d_d);

    if (neighflag != FULL) {
      k_o.template modify<DeviceType>();
      comm->reverse_comm(this); //Coll_vector(q);
      k_o.template sync<DeviceType>();
    }

    // tmp = parallel_dot(d, q, nn);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagQEqDot2>(0,nn),*this,my_dot);
    MPI_Allreduce(&my_dot.v, &dot_sqr.v, 2, MPI_DOUBLE, MPI_SUM, world);
    tmp = dot_sqr;
    if (!(converged & 1))
      alpha[0] = sig_new.v[0] / tmp.v[0];
    if (!(converged & 2))
      alpha[1] = sig_new.v[1] / tmp.v[1];

    sig_old = sig_new;

    // vector_add(s, alpha, d, nn);
    // vector_add(r, -alpha, q, nn);
    // preconditioning: p[j] = r[j] * Hdia_inv[j];
    // sig_new = parallel_dot(r, p, nn);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagQEqDot3>(0,nn),*this,my_dot);
    MPI_Allreduce(&my_dot.v, &dot_sqr.v, 2, MPI_DOUBLE, MPI_SUM, world);
    sig_new = dot_sqr;

    if (!(converged & 1))
      beta[0] = sig_new.v[0] / sig_old.v[0];
    if (!(converged & 2))
      beta[1] = sig_new.v[1] / sig_old.v[1];

    // vector_sum(d, 1., p, beta, d, nn);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqSum1>(0,nn),*this);
  }

  if (loop >= imax && comm->me == 0)
    error->warning(FLERR,fmt::format("Fix qeq/reaxff/kk cg_solve convergence "
                                     "failed after {} iterations at step {}: "
                                     "({}, {})", loop, update->ntimestep,
                                     (sqrt(sig_new.v[0])/b_norm.v[0]), (sqrt(sig_new.v[1])/b_norm.v[1])));

  return loop;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::calculate_q()
{
  F_FLOAT2 sum, sum_all;

  // st_sum = parallel_vector_acc(st, nn);
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagQEqSum2>(0,nn),*this,sum);
  MPI_Allreduce(&sum.v, &sum_all.v, 2, MPI_DOUBLE, MPI_SUM, world);
  const F_FLOAT s_sum = sum_all.v[0];
  const F_FLOAT t_sum = sum_all.v[1];

  // u = s_sum / t_sum;
  delta = s_sum/t_sum;

  // q[i] = s[i] - u * t[i];
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqCalculateQ>(0,nn),*this);
  atomKK->modified(execution_space,Q_MASK);

  pack_flag = 3;
  //comm->forward_comm(this); //Dist_vector(atom->q);
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sparse_matvec_kokkos(typename AT::t_ffloat2_1d &d_xx_in)
{
  d_xx = d_xx_in;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqSparseMatvec1>(0,nn),*this);

  int teamsize;
  int vectorsize;
  int leaguesize;
  if (execution_space == Host) {
    teamsize = 1;
    vectorsize = 1;
    leaguesize = nn;
  } else {
    teamsize = FixQEqReaxFFKokkos<DeviceType>::spmv_teamsize;
    vectorsize = FixQEqReaxFFKokkos<DeviceType>::vectorsize;
    leaguesize = (nn + teamsize - 1) / (teamsize);
  }

  if (neighflag != FULL) {
    int nall = nlocal + atomKK->nghost;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqZeroQGhosts>(atom->nlocal,nall),*this);

    if (need_dup)
      dup_o.reset_except(d_o);

    if (neighflag == HALF)
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceType, TagQEqSparseMatvec2_Half<HALF>>(leaguesize, teamsize, vectorsize), *this);
    else if (neighflag == HALFTHREAD)
      Kokkos::parallel_for(Kokkos::TeamPolicy<DeviceType, TagQEqSparseMatvec2_Half<HALFTHREAD>>(leaguesize, teamsize, vectorsize), *this);

    if (need_dup)
      Kokkos::Experimental::contribute(d_o, dup_o);
  } else // FULL
    Kokkos::parallel_for(Kokkos::TeamPolicy <DeviceType, TagQEqSparseMatvec2_Full>(leaguesize, teamsize, vectorsize), *this);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqSparseMatvec1, const int &ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);
  const auto params_eta = params(itype).eta;
  if (mask[i] & groupbit) {
    if (!(converged & 1))
      d_o(i,0) = params_eta * d_xx(i,0);
    if (!(converged & 2))
      d_o(i,1) = params_eta * d_xx(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqZeroQGhosts, const int &i) const
{
  if (mask[i] & groupbit) {
    if (!(converged & 1))
      d_o(i,0) = 0.0;
    if (!(converged & 2))
      d_o(i,1) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqSparseMatvec2_Half<NEIGHFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagQEqSparseMatvec2_Half<NEIGHFLAG>>::member_type &team) const
{
  int k = team.league_rank() * team.team_size() + team.team_rank();
  if (k < nn) {
    // The q array is duplicated for OpenMP, atomic for GPU, and neither for Serial
    auto v_o = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_o),decltype(ndup_o)>::get(dup_o,ndup_o);
    auto a_o = v_o.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    const int i = d_ilist[k];
    if (mask[i] & groupbit) {
      F_FLOAT2 doitmp;
      const double d_xx_i0 = d_xx(i,0);
      const double d_xx_i1 = d_xx(i,1);

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT2& doi) {
        const int j = d_jlist(jj);
        const auto d_val_jj = d_val(jj);
        if (!(converged & 1)) {
          doi.v[0] += d_val_jj * d_xx(j,0);
          a_o(j,0) += d_val_jj * d_xx_i0;
        }
        if (!(converged & 2)) {
          doi.v[1] += d_val_jj * d_xx(j,1);
          a_o(j,1) += d_val_jj * d_xx_i1;
        }
      }, doitmp);
      Kokkos::single(Kokkos::PerThread(team), [&] () {
        if (!(converged & 1))
          a_o(i,0) += doitmp.v[0];
        if (!(converged & 2))
          a_o(i,1) += doitmp.v[1];
      });
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqSparseMatvec2_Full, const membertype_vec &team) const
{
  int k = team.league_rank() * team.team_size() + team.team_rank();
  if (k < nn) {
    const int i = d_ilist[k];
    if (mask[i] & groupbit) {
      F_FLOAT2 doitmp;
      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT2& doi) {
        const int j = d_jlist(jj);
        const auto d_val_jj = d_val(jj);
        if (!(converged & 1))
          doi.v[0] += d_val_jj * d_xx(j,0);
        if (!(converged & 2))
          doi.v[1] += d_val_jj * d_xx(j,1);
      }, doitmp);
      Kokkos::single(Kokkos::PerThread(team), [&] () {
        if (!(converged & 1))
          d_o(i,0) += doitmp.v[0];
        if (!(converged & 2))
          d_o(i,1) += doitmp.v[1];
      });
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqNorm1, const int &ii, F_FLOAT2& out) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    const auto d_Hdia_inv_i = d_Hdia_inv[i];
    if (!(converged & 1)) {
      d_r(i,0) = 1.0*d_b_st(i,0) + -1.0*d_o(i,0);
      d_d(i,0) = d_r(i,0) * d_Hdia_inv_i;
      out.v[0] += d_b_st(i,0) * d_b_st(i,0);
    }

    if (!(converged & 2)) {
      d_r(i,1) = 1.0*d_b_st(i,1) + -1.0*d_o(i,1);
      d_d(i,1) = d_r(i,1) * d_Hdia_inv_i;
      out.v[1] += d_b_st(i,1) * d_b_st(i,1);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqDot1, const int &ii, F_FLOAT2& out) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    if (!(converged & 1))
      out.v[0] += d_r(i,0) * d_d(i,0);
    if (!(converged & 2))
      out.v[1] += d_r(i,1) * d_d(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqDot2, const int &ii, F_FLOAT2& out) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    if (!(converged & 1))
      out.v[0] += d_d(i,0) * d_o(i,0);
    if (!(converged & 2))
      out.v[1] += d_d(i,1) * d_o(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqDot3, const int &ii, F_FLOAT2& out) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    const auto d_Hdia_inv_i = d_Hdia_inv[i];
    if (!(converged & 1)) {
      const auto alpha_0 = alpha[0];
      d_st(i,0) += alpha_0 * d_d(i,0);
      d_r(i,0) += -alpha_0 * d_o(i,0);
      d_p(i,0) = d_r(i,0) * d_Hdia_inv_i;
      out.v[0] += d_r(i,0) * d_p(i,0);
    }
    if (!(converged & 2)) {
      const auto alpha_1 = alpha[1];
      d_st(i,1) += alpha_1 * d_d(i,1);
      d_r(i,1) += -alpha_1 * d_o(i,1);
      d_p(i,1) = d_r(i,1) * d_Hdia_inv_i;
      out.v[1] += d_r(i,1) * d_p(i,1);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqSum1, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    if (!(converged & 1))
      d_d(i,0) = 1.0 * d_p(i,0) + beta[0] * d_d(i,0);
    if (!(converged & 2))
      d_d(i,1) = 1.0 * d_p(i,1) + beta[1] * d_d(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqSum2, const int &ii, F_FLOAT2& out) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    out.v[0] += d_st(i,0);
    out.v[1] += d_st(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqCalculateQ, const int &ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    q(i) = d_st(i,0) - delta * d_st(i,1);

    for (int k = nprev-1; k > 0; --k) {
      d_s_hist(i,k) = d_s_hist(i,k-1);
      d_t_hist(i,k) = d_t_hist(i,k-1);
    }
    d_s_hist(i,0) = d_st(i,0);
    d_t_hist(i,0) = d_st(i,1);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                                                        int iswap_in, DAT::tdual_xfloat_1d &k_buf,
                                                        int /*pbc_flag*/, int * /*pbc*/)
{
  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  d_buf = k_buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqPackForwardComm>(0,n),*this);
  if (pack_flag == 3) return n;
  else return n*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);

  if (pack_flag == 1) {
    if (!(converged & 1))
      d_buf[i*2] = d_d(j,0);
    if (!(converged & 2))
      d_buf[i*2+1] = d_d(j,1);
  } else if (pack_flag == 2) {
    d_buf[i*2] = d_st(j,0);
    d_buf[i*2+1] = d_st(j,1);
  } else if (pack_flag == 3)
    d_buf[i] = q[j];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  d_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagQEqUnpackForwardComm>(0,n),*this);

  if (pack_flag == 3)
    atomKK->modified(execution_space,Q_MASK); // needed for auto_sync
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqUnpackForwardComm, const int &i) const {
  if (pack_flag == 1) {
    if (!(converged & 1))
      d_d(i+first,0) = d_buf[i*2];
    if (!(converged & 2))
      d_d(i+first,1) = d_buf[i*2+1];
  } else if (pack_flag == 2) {
    d_st(i+first,0) = d_buf[i*2];
    d_st(i+first,1) = d_buf[i*2+1];
  } else if (pack_flag == 3)
    q[i + first] = d_buf[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                    int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  k_d.sync_host();
  if (pack_flag == 1) {
    k_d.sync_host();
    for (m = 0; m < n; m++) {
      if (!(converged & 1))
        buf[m*2] = h_d(list[m],0);
      if (!(converged & 2))
        buf[m*2+1] = h_d(list[m],1);
    }
  } else if (pack_flag == 2) {
    k_st.sync_host();
    for (m = 0; m < n; m++) {
      buf[m*2] = h_st(list[m],0);
      buf[m*2+1] = h_st(list[m],1);
    }
  } else if (pack_flag == 3) {
    atomKK->sync(Host,Q_MASK);
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  }

  if (pack_flag == 3) return n;
  else return n*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1) {
    k_d.sync_host();
    for (m = 0, i = first; m < n; m++, i++) {
      if (!(converged & 1))
        h_d(i,0) = buf[m*2];
      if (!(converged & 2))
        h_d(i,1) = buf[m*2+1];
    }
    k_d.modify_host();
  } else if (pack_flag == 2) {
    k_st.sync_host();
    for (m = 0, i = first; m < n; m++, i++) {
      h_st(i,0) = buf[m*2];
      h_st(i,1) = buf[m*2+1];
    }
    k_st.modify_host();
  } else if (pack_flag == 3) {
    atomKK->sync(Host,Q_MASK);
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
    atomKK->modified(Host,Q_MASK);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  k_o.sync_host();
  for (m = 0, i = first; m < n; m++, i++) {
    if (!(converged & 1))
      buf[m*2] = h_o(i,0);
    if (!(converged & 2))
      buf[m*2+1] = h_o(i,1);
  }
  return n*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  k_o.sync_host();
  for (int m = 0; m < n; m++) {
    if (!(converged & 1))
      h_o(list[m],0) += buf[m*2];
    if (!(converged & 2))
      h_o(list[m],1) += buf[m*2+1];
  }
  k_o.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixQEqReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes;

  bytes = atom->nmax*nprev*2 * sizeof(F_FLOAT); // s_hist & t_hist
  bytes += (double)atom->nmax*8 * sizeof(F_FLOAT); // storage
  bytes += (double)n_cap*2 * sizeof(int); // matrix...
  bytes += (double)m_cap * sizeof(int);
  bytes += (double)m_cap * sizeof(F_FLOAT);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::grow_arrays(int nmax)
{
  k_s_hist.template sync<LMPHostType>();
  k_t_hist.template sync<LMPHostType>();

  k_s_hist.template modify<LMPHostType>(); // force reallocation on host
  k_t_hist.template modify<LMPHostType>();

  memoryKK->grow_kokkos(k_s_hist,s_hist,nmax,nprev,"qeq:s_hist");
  memoryKK->grow_kokkos(k_t_hist,t_hist,nmax,nprev,"qeq:t_hist");

  d_s_hist = k_s_hist.template view<DeviceType>();
  d_t_hist = k_t_hist.template view<DeviceType>();

  k_s_hist.template modify<LMPHostType>();
  k_t_hist.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  k_s_hist.template sync<LMPHostType>();
  k_t_hist.template sync<LMPHostType>();

  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }

  k_s_hist.template modify<LMPHostType>();
  k_t_hist.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   sort local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  // always sort on the device

  k_s_hist.sync_device();
  k_t_hist.sync_device();

  Sorter.sort(LMPDeviceType(), k_s_hist.d_view);
  Sorter.sort(LMPDeviceType(), k_t_hist.d_view);

  k_s_hist.modify_device();
  k_t_hist.modify_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqPackExchange, const int &mysend) const {
  const int i = d_exchange_sendlist(mysend);

  for (int m = 0; m < nprev; m++) d_buf(mysend*nprev*2 + m) = d_s_hist(i,m);
  for (int m = 0; m < nprev; m++) d_buf(mysend*nprev*2 + nprev+m) = d_t_hist(i,m);

  const int j = d_copylist(mysend);

  if (j > -1) {
    for (int m = 0; m < nprev; m++) d_s_hist(i,m) = d_s_hist(j,m);
    for (int m = 0; m < nprev; m++) d_t_hist(i,m) = d_t_hist(j,m);
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

  k_s_hist.template sync<DeviceType>();
  k_t_hist.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqPackExchange>(0,nsend),*this);

  copymode = 0;

  k_s_hist.template modify<DeviceType>();
  k_t_hist.template modify<DeviceType>();

  return nsend*nprev*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxFFKokkos<DeviceType>::operator()(TagQEqUnpackExchange, const int &i) const
{
  int index = d_indices(i);

  if (index > -1) {
    for (int m = 0; m < nprev; m++) d_s_hist(index,m) = d_buf(i*nprev*2 + m);
    for (int m = 0; m < nprev; m++) d_t_hist(index,m) = d_buf(i*nprev*2 + nprev+m);
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

  k_s_hist.template sync<DeviceType>();
  k_t_hist.template sync<DeviceType>();

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagQEqUnpackExchange>(0,nrecv),*this);

  copymode = 0;

  k_s_hist.template modify<DeviceType>();
  k_t_hist.template modify<DeviceType>();
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_s_hist.template sync<LMPHostType>();
  k_t_hist.template sync<LMPHostType>();

  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];

  k_s_hist.template modify<LMPHostType>();
  k_t_hist.template modify<LMPHostType>();

  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxFFKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  k_s_hist.template sync<LMPHostType>();
  k_t_hist.template sync<LMPHostType>();

  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];

  k_s_hist.template modify<LMPHostType>();
  k_t_hist.template modify<LMPHostType>();

  return nprev*2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxFFKokkos<DeviceType>::get_chi_field()
{
  atomKK->sync(Host,X_MASK|MASK_MASK|IMAGE_MASK);
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
