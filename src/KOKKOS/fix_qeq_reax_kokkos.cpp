/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ray Shan (SNL), Stan Moore (SNL),
                          Kamesh Arumugam (NVIDIA)
------------------------------------------------------------------------- */

#include "fix_qeq_reax_kokkos.h"
#include <cmath>
#include "kokkos.h"
#include "atom.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "memory_kokkos.h"
#include "error.h"
#include "pair_reaxc_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 0.0001
#define EV_TO_KCAL_PER_MOL 14.4

#define TEAMSIZE 128

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxKokkos<DeviceType>::
FixQEqReaxKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixQEqReax(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | Q_MASK | TYPE_MASK | TAG_MASK;
  datamask_modify = Q_MASK | X_MASK;

  nmax = nmax = m_cap = 0;
  allocated_flag = 0;
  nprev = 4;

  memory->destroy(s_hist);
  memory->destroy(t_hist);
  grow_arrays(atom->nmax);

  d_mfill_offset = typename AT::t_int_scalar("qeq/kk:mfill_offset");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixQEqReaxKokkos<DeviceType>::~FixQEqReaxKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_s_hist,s_hist);
  memoryKK->destroy_kokkos(k_t_hist,t_hist);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::init()
{
  atomKK->k_q.modify<LMPHostType>();
  atomKK->k_q.sync<DeviceType>();

  FixQEqReax::init();

  neighflag = lmp->kokkos->neighflag_qeq;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL) {
    neighbor->requests[irequest]->fix = 1;
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else { //if (neighflag == HALF || neighflag == HALFTHREAD)
    neighbor->requests[irequest]->fix = 1;
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
    neighbor->requests[irequest]->ghost = 1;
  }

  int ntypes = atom->ntypes;
  k_params = Kokkos::DualView<params_qeq*,Kokkos::LayoutRight,DeviceType>
    ("FixQEqReax::params",ntypes+1);
  params = k_params.template view<DeviceType>();

  for (n = 1; n <= ntypes; n++) {
    k_params.h_view(n).chi = chi[n];
    k_params.h_view(n).eta = eta[n];
    k_params.h_view(n).gamma = gamma[n];
  }
  k_params.template modify<LMPHostType>();

  cutsq = swb * swb;

  init_shielding_k();
  init_hist();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::init_shielding_k()
{
  int i,j;
  int ntypes = atom->ntypes;

  k_shield = DAT::tdual_ffloat_2d("qeq/kk:shield",ntypes+1,ntypes+1);
  d_shield = k_shield.template view<DeviceType>();

  for( i = 1; i <= ntypes; ++i )
    for( j = 1; j <= ntypes; ++j )
      k_shield.h_view(i,j) = pow( gamma[i] * gamma[j], -1.5 );

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
void FixQEqReaxKokkos<DeviceType>::init_hist()
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
void FixQEqReaxKokkos<DeviceType>::setup_pre_force(int vflag)
{
  //neighbor->build_one(list);

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::pre_force(int vflag)
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
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;

  k_params.template sync<DeviceType>();
  k_shield.template sync<DeviceType>();
  k_tap.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  copymode = 1;

  // allocate

  allocate_array();

  // get max number of neighbor

  if (!allocated_flag || update->ntimestep == neighbor->lastcall)
    allocate_matrix();

  // compute_H

  if (lmp->kokkos->ngpus == 0) { // CPU
    if (neighflag == FULL) {
      FixQEqReaxKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(this);
      Kokkos::parallel_scan(inum,computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixQEqReaxKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(this);
      Kokkos::parallel_scan(inum,computeH_functor);
    }
  } else { // GPU, use teams
    Kokkos::deep_copy(d_mfill_offset,0);

    int vector_length = 32;
    int atoms_per_team = 4;
    int num_teams = inum / atoms_per_team + (inum % atoms_per_team ? 1 : 0);

    Kokkos::TeamPolicy<DeviceType> policy(num_teams, atoms_per_team,
                                          vector_length);
    if (neighflag == FULL) {
      FixQEqReaxKokkosComputeHFunctor<DeviceType, FULL> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    } else { // HALF and HALFTHREAD are the same
      FixQEqReaxKokkosComputeHFunctor<DeviceType, HALF> computeH_functor(
          this, atoms_per_team, vector_length);
      Kokkos::parallel_for(policy, computeH_functor);
    }
  }

  // init_matvec

  k_s_hist.template sync<DeviceType>();
  k_t_hist.template sync<DeviceType>();
  FixQEqReaxKokkosMatVecFunctor<DeviceType> matvec_functor(this);
  Kokkos::parallel_for(inum,matvec_functor);

  // comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 2;
  k_s.template modify<DeviceType>();
  k_s.template sync<LMPHostType>();
  comm->forward_comm_fix(this);
  k_s.template modify<LMPHostType>();
  k_s.template sync<DeviceType>();

  // comm->forward_comm_fix(this); //Dist_vector( t );
  pack_flag = 3;
  k_t.template modify<DeviceType>();
  k_t.template sync<LMPHostType>();
  comm->forward_comm_fix(this);
  k_t.template modify<LMPHostType>();
  k_t.template sync<DeviceType>();

  need_dup = lmp->kokkos->need_dup<DeviceType>();

  if (need_dup)
    dup_o = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated> (d_o); // allocate duplicated memory
  else
    ndup_o = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated> (d_o);

  // 1st cg solve over b_s, s

  cg_solve1();

  // 2nd cg solve over b_t, t

  cg_solve2();

  // calculate_Q();

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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::num_neigh_item(int ii, int &maxneigh) const
{
  const int i = d_ilist[ii];
  maxneigh += d_numneigh[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::allocate_matrix()
{
  const int inum = list->inum;

  nmax = atom->nmax;

  // determine the total space for the H matrix

  m_cap = 0;
  FixQEqReaxKokkosNumNeighFunctor<DeviceType> neigh_functor(this);
  Kokkos::parallel_reduce(inum,neigh_functor,m_cap);

  d_firstnbr = typename AT::t_int_1d("qeq/kk:firstnbr",nmax);
  d_numnbrs = typename AT::t_int_1d("qeq/kk:numnbrs",nmax);
  d_jlist = typename AT::t_int_1d("qeq/kk:jlist",m_cap);
  d_val = typename AT::t_ffloat_1d("qeq/kk:val",m_cap);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::allocate_array()
{
  if (atom->nmax > nmax) {
    nmax = atom->nmax;

    k_o = DAT::tdual_ffloat_1d("qeq/kk:h_o",nmax);
    d_o = k_o.template view<DeviceType>();
    h_o = k_o.h_view;

    d_Hdia_inv = typename AT::t_ffloat_1d("qeq/kk:h_Hdia_inv",nmax);

    d_b_s = typename AT::t_ffloat_1d("qeq/kk:h_b_s",nmax);

    d_b_t = typename AT::t_ffloat_1d("qeq/kk:h_b_t",nmax);

    k_s = DAT::tdual_ffloat_1d("qeq/kk:h_s",nmax);
    d_s = k_s.template view<DeviceType>();
    h_s = k_s.h_view;

    k_t = DAT::tdual_ffloat_1d("qeq/kk:h_t",nmax);
    d_t = k_t.template view<DeviceType>();
    h_t = k_t.h_view;

    d_p = typename AT::t_ffloat_1d("qeq/kk:h_p",nmax);

    d_r = typename AT::t_ffloat_1d("qeq/kk:h_r",nmax);

    k_d = DAT::tdual_ffloat_1d("qeq/kk:h_d",nmax);
    d_d = k_d.template view<DeviceType>();
    h_d = k_d.h_view;
  }

  // init_storage
  const int ignum = atom->nlocal + atom->nghost;
  FixQEqReaxKokkosZeroFunctor<DeviceType> zero_functor(this);
  Kokkos::parallel_for(ignum,zero_functor);

}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::zero_item(int ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_s[i] = -params(itype).chi;
    d_b_t[i] = -1.0;
    d_s[i] = 0.0;
    d_t[i] = 0.0;
    d_p[i] = 0.0;
    d_o[i] = 0.0;
    d_r[i] = 0.0;
    d_d[i] = 0.0;
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::compute_h_item(int ii, int &m_fill, const bool &final) const
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
void FixQEqReaxKokkos<DeviceType>::compute_h_team(
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
      (firstatom + atoms_per_team < inum) ? (firstatom + atoms_per_team) : inum;

  // kokkos-thread-0 is used to load info from global memory into scratch space
  if (team.team_rank() == 0) {

    // copy atom indices from d_ilist[firstatom:lastatom] to scratch space s_ilist[0:atoms_per_team]
    // copy # of neighbor atoms for all the atoms with indices in d_ilist[firstatom:lastatom] from d_numneigh to scratch space s_numneigh[0:atoms_per_team]
    // calculate total number of neighbor atoms for all atoms assigned to the current team of threads (Note - Total # of neighbor atoms here provides the
    // upper bound space requirement to store the H matrix values corresponding to the atoms with indices in d_ilist[firstatom:lastatom])

    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team, atoms_per_team),
                          [&](const int &idx, int &totalnbrs, bool final) {
                            int ii = firstatom + idx;

                            if (ii < inum) {
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

        if (ii < inum) {
          const int i = s_ilist[idx];

          if (mask[i] & groupbit) {
            const X_FLOAT xtmp = x(i, 0);
            const X_FLOAT ytmp = x(i, 1);
            const X_FLOAT ztmp = x(i, 2);
            const int itype = type(i);
            const tagint itag = tag(i);
            const int jnum = s_numnbrs[idx];

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
double FixQEqReaxKokkos<DeviceType>::calculate_H_k(const F_FLOAT &r, const F_FLOAT &shld) const
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
  denom = pow(denom,0.3333333333333);

  return taper * EV_TO_KCAL_PER_MOL / denom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::matvec_item(int ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);

  if (mask[i] & groupbit) {
    d_Hdia_inv[i] = 1.0 / params(itype).eta;
    d_b_s[i] = -params(itype).chi;
    d_b_t[i] = -1.0;
    d_t[i] = d_t_hist(i,2) + 3*(d_t_hist(i,0) - d_t_hist(i,1));
    d_s[i] = 4*(d_s_hist(i,0)+d_s_hist(i,2))-(6*d_s_hist(i,1)+d_s_hist(i,3));
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::cg_solve1()
// b = b_s, x = s;
{
  const int inum = list->inum;
  F_FLOAT tmp, sig_old, b_norm;

  const int teamsize = TEAMSIZE;

  // sparse_matvec( &H, x, q );
  FixQEqReaxKokkosSparse12Functor<DeviceType> sparse12_functor(this);
  Kokkos::parallel_for(inum,sparse12_functor);
  if (neighflag != FULL) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagZeroQGhosts>(nlocal,nlocal+atom->nghost),*this);
    if (neighflag == HALF) {
      FixQEqReaxKokkosSparse13Functor<DeviceType,HALF> sparse13_functor(this);
      Kokkos::parallel_for(inum,sparse13_functor);
    } else if (neighflag == HALFTHREAD) {
      FixQEqReaxKokkosSparse13Functor<DeviceType,HALFTHREAD> sparse13_functor(this);
      Kokkos::parallel_for(inum,sparse13_functor);
    }
    if (need_dup)
      Kokkos::Experimental::contribute(d_o, dup_o);
  } else {
    Kokkos::parallel_for(Kokkos::TeamPolicy <DeviceType, TagSparseMatvec1> (inum, teamsize), *this);
  }

  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    k_o.template sync<LMPHostType>();
    comm->reverse_comm_fix(this); //Coll_vector( q );
    k_o.template modify<LMPHostType>();
    k_o.template sync<DeviceType>();
  }

  // vector_sum( r , 1.,  b, -1., q, nn );
  // preconditioning: d[j] = r[j] * Hdia_inv[j];
  // b_norm = parallel_norm( b, nn );
  F_FLOAT my_norm = 0.0;
  FixQEqReaxKokkosNorm1Functor<DeviceType> norm1_functor(this);
  Kokkos::parallel_reduce(inum,norm1_functor,my_norm);
  F_FLOAT norm_sqr = 0.0;
  MPI_Allreduce( &my_norm, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  b_norm = sqrt(norm_sqr);

  // sig_new = parallel_dot( r, d, nn);
  F_FLOAT my_dot = 0.0;
  FixQEqReaxKokkosDot1Functor<DeviceType> dot1_functor(this);
  Kokkos::parallel_reduce(inum,dot1_functor,my_dot);
  F_FLOAT dot_sqr = 0.0;
  MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  F_FLOAT sig_new = dot_sqr;

  int loop;
  const int loopmax = 200;
  for (loop = 1; loop < loopmax & sqrt(sig_new)/b_norm > tolerance; loop++) {

    // comm->forward_comm_fix(this); //Dist_vector( d );
    pack_flag = 1;
    k_d.template modify<DeviceType>();
    k_d.template sync<LMPHostType>();
    comm->forward_comm_fix(this);
    k_d.template modify<LMPHostType>();
    k_d.template sync<DeviceType>();

    // sparse_matvec( &H, d, q );
    FixQEqReaxKokkosSparse22Functor<DeviceType> sparse22_functor(this);
    Kokkos::parallel_for(inum,sparse22_functor);
    if (neighflag != FULL) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagZeroQGhosts>(nlocal,nlocal+atom->nghost),*this);
      if (need_dup)
        dup_o.reset_except(d_o);
      if (neighflag == HALF) {
        FixQEqReaxKokkosSparse23Functor<DeviceType,HALF> sparse23_functor(this);
        Kokkos::parallel_for(inum,sparse23_functor);
      } else if (neighflag == HALFTHREAD) {
        FixQEqReaxKokkosSparse23Functor<DeviceType,HALFTHREAD> sparse23_functor(this);
        Kokkos::parallel_for(inum,sparse23_functor);
      }
      if (need_dup)
        Kokkos::Experimental::contribute(d_o, dup_o);
    } else {
      Kokkos::parallel_for(Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2> (inum, teamsize), *this);
    }

    if (neighflag != FULL) {
      k_o.template modify<DeviceType>();
      k_o.template sync<LMPHostType>();
      comm->reverse_comm_fix(this); //Coll_vector( q );
      k_o.template modify<LMPHostType>();
      k_o.template sync<DeviceType>();
    }

    // tmp = parallel_dot( d, q, nn);
    my_dot = dot_sqr = 0.0;
    FixQEqReaxKokkosDot2Functor<DeviceType> dot2_functor(this);
    Kokkos::parallel_reduce(inum,dot2_functor,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    tmp = dot_sqr;

    alpha = sig_new / tmp;

    sig_old = sig_new;

    // vector_add( s, alpha, d, nn );
    // vector_add( r, -alpha, q, nn );
    my_dot = dot_sqr = 0.0;
    FixQEqReaxKokkosPrecon1Functor<DeviceType> precon1_functor(this);
    Kokkos::parallel_for(inum,precon1_functor);
    // preconditioning: p[j] = r[j] * Hdia_inv[j];
    // sig_new = parallel_dot( r, p, nn);
    FixQEqReaxKokkosPreconFunctor<DeviceType> precon_functor(this);
    Kokkos::parallel_reduce(inum,precon_functor,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    sig_new = dot_sqr;

    beta = sig_new / sig_old;

    // vector_sum( d, 1., p, beta, d, nn );
    FixQEqReaxKokkosVecSum2Functor<DeviceType> vecsum2_functor(this);
    Kokkos::parallel_for(inum,vecsum2_functor);
  }

  if (loop >= loopmax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq/reax cg_solve1 convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step: %f",loop,update->ntimestep,sqrt(sig_new)/b_norm);
    error->warning(FLERR,str);
    //error->all(FLERR,str);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::cg_solve2()
// b = b_t, x = t;
{
  const int inum = list->inum;
  F_FLOAT tmp, sig_old, b_norm;

  const int teamsize = TEAMSIZE;

  // sparse_matvec( &H, x, q );
  FixQEqReaxKokkosSparse32Functor<DeviceType> sparse32_functor(this);
  Kokkos::parallel_for(inum,sparse32_functor);
  if (neighflag != FULL) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagZeroQGhosts>(nlocal,nlocal+atom->nghost),*this);
    if (need_dup)
      dup_o.reset_except(d_o);
    if (neighflag == HALF) {
      FixQEqReaxKokkosSparse33Functor<DeviceType,HALF> sparse33_functor(this);
      Kokkos::parallel_for(inum,sparse33_functor);
    } else if (neighflag == HALFTHREAD) {
      FixQEqReaxKokkosSparse33Functor<DeviceType,HALFTHREAD> sparse33_functor(this);
      Kokkos::parallel_for(inum,sparse33_functor);
    }
    if (need_dup)
      Kokkos::Experimental::contribute(d_o, dup_o);
  } else {
    Kokkos::parallel_for(Kokkos::TeamPolicy <DeviceType, TagSparseMatvec3> (inum, teamsize), *this);
  }

  if (neighflag != FULL) {
    k_o.template modify<DeviceType>();
    k_o.template sync<LMPHostType>();
    comm->reverse_comm_fix(this); //Coll_vector( q );
    k_o.template modify<LMPHostType>();
    k_o.template sync<DeviceType>();
  }

  // vector_sum( r , 1.,  b, -1., q, nn );
  // preconditioning: d[j] = r[j] * Hdia_inv[j];
  // b_norm = parallel_norm( b, nn );
  F_FLOAT my_norm = 0.0;
  FixQEqReaxKokkosNorm2Functor<DeviceType> norm2_functor(this);
  Kokkos::parallel_reduce(inum,norm2_functor,my_norm);
  F_FLOAT norm_sqr = 0.0;
  MPI_Allreduce( &my_norm, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  b_norm = sqrt(norm_sqr);

  // sig_new = parallel_dot( r, d, nn);
  F_FLOAT my_dot = 0.0;
  FixQEqReaxKokkosDot1Functor<DeviceType> dot1_functor(this);
  Kokkos::parallel_reduce(inum,dot1_functor,my_dot);
  F_FLOAT dot_sqr = 0.0;
  MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
  F_FLOAT sig_new = dot_sqr;

  int loop;
  const int loopmax = 200;
  for (loop = 1; loop < loopmax & sqrt(sig_new)/b_norm > tolerance; loop++) {

    // comm->forward_comm_fix(this); //Dist_vector( d );
    pack_flag = 1;
    k_d.template modify<DeviceType>();
    k_d.template sync<LMPHostType>();
    comm->forward_comm_fix(this);
    k_d.template modify<LMPHostType>();
    k_d.template sync<DeviceType>();

    // sparse_matvec( &H, d, q );
    FixQEqReaxKokkosSparse22Functor<DeviceType> sparse22_functor(this);
    Kokkos::parallel_for(inum,sparse22_functor);
    if (neighflag != FULL) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagZeroQGhosts>(nlocal,nlocal+atom->nghost),*this);
      if (need_dup)
        dup_o.reset_except(d_o);
      if (neighflag == HALF) {
        FixQEqReaxKokkosSparse23Functor<DeviceType,HALF> sparse23_functor(this);
        Kokkos::parallel_for(inum,sparse23_functor);
      } else if (neighflag == HALFTHREAD) {
        FixQEqReaxKokkosSparse23Functor<DeviceType,HALFTHREAD> sparse23_functor(this);
        Kokkos::parallel_for(inum,sparse23_functor);
      }
      if (need_dup)
        Kokkos::Experimental::contribute(d_o, dup_o);
    } else {
      Kokkos::parallel_for(Kokkos::TeamPolicy <DeviceType, TagSparseMatvec2> (inum, teamsize), *this);
    }

    if (neighflag != FULL) {
      k_o.template modify<DeviceType>();
      k_o.template sync<LMPHostType>();
      comm->reverse_comm_fix(this); //Coll_vector( q );
      k_o.template modify<LMPHostType>();
      k_o.template sync<DeviceType>();
    }

    // tmp = parallel_dot( d, q, nn);
    my_dot = dot_sqr = 0.0;
    FixQEqReaxKokkosDot2Functor<DeviceType> dot2_functor(this);
    Kokkos::parallel_reduce(inum,dot2_functor,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    tmp = dot_sqr;

    alpha = sig_new / tmp;

    sig_old = sig_new;

    // vector_add( t, alpha, d, nn );
    // vector_add( r, -alpha, q, nn );
    my_dot = dot_sqr = 0.0;
    FixQEqReaxKokkosPrecon2Functor<DeviceType> precon2_functor(this);
    Kokkos::parallel_for(inum,precon2_functor);
    // preconditioning: p[j] = r[j] * Hdia_inv[j];
    // sig_new = parallel_dot( r, p, nn);
    FixQEqReaxKokkosPreconFunctor<DeviceType> precon_functor(this);
    Kokkos::parallel_reduce(inum,precon_functor,my_dot);
    MPI_Allreduce( &my_dot, &dot_sqr, 1, MPI_DOUBLE, MPI_SUM, world );
    sig_new = dot_sqr;

    beta = sig_new / sig_old;

    // vector_sum( d, 1., p, beta, d, nn );
    FixQEqReaxKokkosVecSum2Functor<DeviceType> vecsum2_functor(this);
    Kokkos::parallel_for(inum,vecsum2_functor);
  }

  if (loop >= loopmax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix qeq/reax cg_solve2 convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step: %f",loop,update->ntimestep,sqrt(sig_new)/b_norm);
    error->warning(FLERR,str);
    //error->all(FLERR,str);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::calculate_q()
{
  F_FLOAT sum, sum_all;
  const int inum = list->inum;

  // s_sum = parallel_vector_acc( s, nn );
  sum = sum_all = 0.0;
  FixQEqReaxKokkosVecAcc1Functor<DeviceType> vecacc1_functor(this);
  Kokkos::parallel_reduce(inum,vecacc1_functor,sum);
  MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, world );
  const F_FLOAT s_sum = sum_all;

  // t_sum = parallel_vector_acc( t, nn);
  sum = sum_all = 0.0;
  FixQEqReaxKokkosVecAcc2Functor<DeviceType> vecacc2_functor(this);
  Kokkos::parallel_reduce(inum,vecacc2_functor,sum);
  MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, world );
  const F_FLOAT t_sum = sum_all;

  // u = s_sum / t_sum;
  delta = s_sum/t_sum;

  // q[i] = s[i] - u * t[i];
  FixQEqReaxKokkosCalculateQFunctor<DeviceType> calculateQ_functor(this);
  Kokkos::parallel_for(inum,calculateQ_functor);

  pack_flag = 4;
  //comm->forward_comm_fix( this ); //Dist_vector( atom->q );
  atomKK->k_q.modify<DeviceType>();
  atomKK->k_q.sync<LMPHostType>();
  comm->forward_comm_fix(this);
  atomKK->k_q.modify<LMPHostType>();
  atomKK->k_q.sync<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse12_item(int ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);
  if (mask[i] & groupbit) {
    d_o[i] = params(itype).eta * d_s[i];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse13_item(int ii) const
{
  // The q array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_o = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_o),decltype(ndup_o)>::get(dup_o,ndup_o);
  auto a_o = v_o.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    F_FLOAT tmp = 0.0;
    for(int jj = d_firstnbr[i]; jj < d_firstnbr[i] + d_numnbrs[i]; jj++) {
      const int j = d_jlist(jj);
      tmp += d_val(jj) * d_s[j];
      a_o[j] += d_val(jj) * d_s[i];
    }
    a_o[i] += tmp;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::operator() (TagSparseMatvec1, const membertype1 &team) const
{
  const int i = d_ilist[team.league_rank()];
  if (mask[i] & groupbit) {
    F_FLOAT doitmp;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT &doi) {
      const int j = d_jlist(jj);
      doi += d_val(jj) * d_s[j];
    }, doitmp);
    Kokkos::single(Kokkos::PerTeam(team), [&] () {d_o[i] += doitmp;});
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse22_item(int ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);
  if (mask[i] & groupbit) {
    d_o[i] = params(itype).eta * d_d[i];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse23_item(int ii) const
{
  // The q array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_o = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_o),decltype(ndup_o)>::get(dup_o,ndup_o);
  auto a_o = v_o.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    F_FLOAT tmp = 0.0;
    for(int jj = d_firstnbr[i]; jj < d_firstnbr[i] + d_numnbrs[i]; jj++) {
      const int j = d_jlist(jj);
      tmp += d_val(jj) * d_d[j];
      a_o[j] += d_val(jj) * d_d[i];
    }
    a_o[i] += tmp;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::operator() (TagSparseMatvec2, const membertype2 &team) const
{
  const int i = d_ilist[team.league_rank()];
  if (mask[i] & groupbit) {
    F_FLOAT doitmp;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT &doi) {
      const int j = d_jlist(jj);
      doi += d_val(jj) * d_d[j];
    }, doitmp);
    Kokkos::single(Kokkos::PerTeam(team), [&] () {d_o[i] += doitmp; });
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::operator() (TagZeroQGhosts, const int &i) const
{
  if (mask[i] & groupbit)
    d_o[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse32_item(int ii) const
{
  const int i = d_ilist[ii];
  const int itype = type(i);
  if (mask[i] & groupbit)
    d_o[i] = params(itype).eta * d_t[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::sparse33_item(int ii) const
{
  // The q array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_o = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_o),decltype(ndup_o)>::get(dup_o,ndup_o);
  auto a_o = v_o.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    F_FLOAT tmp = 0.0;
    for(int jj = d_firstnbr[i]; jj < d_firstnbr[i] + d_numnbrs[i]; jj++) {
      const int j = d_jlist(jj);
      tmp += d_val(jj) * d_t[j];
      a_o[j] += d_val(jj) * d_t[i];
    }
    a_o[i] += tmp;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::operator() (TagSparseMatvec3, const membertype3 &team) const
{
  const int i = d_ilist[team.league_rank()];
  if (mask[i] & groupbit) {
    F_FLOAT doitmp;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, d_firstnbr[i], d_firstnbr[i] + d_numnbrs[i]), [&] (const int &jj, F_FLOAT &doi) {
      const int j = d_jlist(jj);
      doi += d_val(jj) * d_t[j];
    }, doitmp);
    Kokkos::single(Kokkos::PerTeam(team), [&] () {d_o[i] += doitmp;});
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::vecsum2_item(int ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit)
    d_d[i] = 1.0 * d_p[i] + beta * d_d[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::norm1_item(int ii) const
{
  F_FLOAT tmp = 0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_r[i] = 1.0*d_b_s[i] + -1.0*d_o[i];
    d_d[i] = d_r[i] * d_Hdia_inv[i];
    tmp = d_b_s[i] * d_b_s[i];
  }
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::norm2_item(int ii) const
{
  F_FLOAT tmp = 0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_r[i] = 1.0*d_b_t[i] + -1.0*d_o[i];
    d_d[i] = d_r[i] * d_Hdia_inv[i];
    tmp = d_b_t[i] * d_b_t[i];
  }
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::dot1_item(int ii) const
{
  F_FLOAT tmp = 0.0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit)
    tmp = d_r[i] * d_d[i];
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::dot2_item(int ii) const
{
  double tmp = 0.0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    tmp = d_d[i] * d_o[i];
  }
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::precon1_item(int ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_s[i] += alpha * d_d[i];
    d_r[i] += -alpha * d_o[i];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::precon2_item(int ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_t[i] += alpha * d_d[i];
    d_r[i] += -alpha * d_o[i];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::precon_item(int ii) const
{
  F_FLOAT tmp = 0.0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    d_p[i] = d_r[i] * d_Hdia_inv[i];
    tmp = d_r[i] * d_p[i];
  }
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::vecacc1_item(int ii) const
{
  F_FLOAT tmp = 0.0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit)
    tmp = d_s[i];
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double FixQEqReaxKokkos<DeviceType>::vecacc2_item(int ii) const
{
  F_FLOAT tmp = 0.0;
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    tmp = d_t[i];
  }
  return tmp;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixQEqReaxKokkos<DeviceType>::calculate_q_item(int ii) const
{
  const int i = d_ilist[ii];
  if (mask[i] & groupbit) {
    q(i) = d_s[i] - delta * d_t[i];

    for (int k = nprev-1; k > 0; --k) {
      d_s_hist(i,k) = d_s_hist(i,k-1);
      d_t_hist(i,k) = d_t_hist(i,k-1);
    }
    d_s_hist(i,0) = d_s[i];
    d_t_hist(i,0) = d_t[i];
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int m;

  if (pack_flag == 1)
    for(m = 0; m < n; m++) buf[m] = h_d[list[m]];
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = h_s[list[m]];
  else if( pack_flag == 3 )
    for(m = 0; m < n; m++) buf[m] = h_t[list[m]];
  else if( pack_flag == 4 )
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];

  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for(m = 0, i = first; m < n; m++, i++) h_d[i] = buf[m];
  else if( pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) h_s[i] = buf[m];
  else if( pack_flag == 3)
    for(m = 0, i = first; m < n; m++, i++) h_t[i] = buf[m];
  else if( pack_flag == 4)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for(m = 0, i = first; m < n; m++, i++) {
    buf[m] = h_o[i];
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  for(int m = 0; m < n; m++) {
    h_o[list[m]] += buf[m];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double FixQEqReaxKokkos<DeviceType>::memory_usage()
{
  double bytes;

  bytes = atom->nmax*nprev*2 * sizeof(F_FLOAT); // s_hist & t_hist
  bytes += atom->nmax*8 * sizeof(F_FLOAT); // storage
  bytes += n_cap*2 * sizeof(int); // matrix...
  bytes += m_cap * sizeof(int);
  bytes += m_cap * sizeof(F_FLOAT);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void FixQEqReaxKokkos<DeviceType>::grow_arrays(int nmax)
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
void FixQEqReaxKokkos<DeviceType>::copy_arrays(int i, int j, int delflag)
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
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_s_hist.template sync<LMPHostType>();
  k_t_hist.template sync<LMPHostType>();

  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

template<class DeviceType>
int FixQEqReaxKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];

  k_s_hist.template modify<LMPHostType>();
  k_t_hist.template modify<LMPHostType>();

  return nprev*2;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixQEqReaxKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixQEqReaxKokkos<LMPHostType>;
#endif
}
