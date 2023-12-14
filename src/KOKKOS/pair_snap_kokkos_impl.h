// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL),
                         Evan Weinberg (NVIDIA)
------------------------------------------------------------------------- */

#include "pair_snap_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"
#include "sna.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#define MAXLINE 1024
#define MAXWORD 3

namespace LAMMPS_NS {

// Outstanding issues with quadratic term
// 1. there seems to a problem with compute_optimized energy calc
// it does not match compute_regular, even when quadratic coeffs = 0

//static double t1 = 0.0;
//static double t2 = 0.0;
//static double t3 = 0.0;
//static double t4 = 0.0;
//static double t5 = 0.0;
//static double t6 = 0.0;
//static double t7 = 0.0;
/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
PairSNAPKokkos<DeviceType, real_type, vector_length>::PairSNAPKokkos(LAMMPS *lmp) : PairSNAP(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
PairSNAPKokkos<DeviceType, real_type, vector_length>::~PairSNAPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::init_style()
{
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR,"Pair style snap/kk can currently only run on a single "
                         "CPU thread");

    PairSNAP::init_style();
    return;
  }

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair snap/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& max_neighs) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }
};

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::compute(int eflag_in, int vflag_in)
{
  if (host_flag) {
    atomKK->sync(Host,X_MASK|F_MASK|TYPE_MASK);
    PairSNAP::compute(eflag_in,vflag_in);
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
  int newton_pair = force->newton_pair;
  if (newton_pair == false)
    error->all(FLERR,"PairSNAPKokkos requires 'newton on'");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
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

  /*
  for (int i = 0; i < nlocal; i++) {
    typename t_neigh_list::t_neighs neighs_i = neigh_list.get_neighs(i);
    const int num_neighs = neighs_i.get_num_neighs();
    if (max_neighs<num_neighs) max_neighs = num_neighs;
  }*/
  max_neighs = 0;
  Kokkos::parallel_reduce("PairSNAPKokkos::find_max_neighs",inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(max_neighs));

  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 32;//max_neighs;

  if (beta_max < inum) {
    beta_max = inum;
    MemKK::realloc_kokkos(d_beta,"PairSNAPKokkos:beta",ncoeff,inum);
    if (!host_flag)
      MemKK::realloc_kokkos(d_beta_pack,"PairSNAPKokkos:beta_pack",vector_length,ncoeff,(inum + vector_length - 1) / vector_length);
    MemKK::realloc_kokkos(d_ninside,"PairSNAPKokkos:ninside",inum);
  }

  chunk_size = MIN(chunksize,inum); // "chunksize" variable is set by user
  chunk_offset = 0;

  snaKK.grow_rij(chunk_size,max_neighs);

  EV_FLOAT ev;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    if (host_flag)
    {
      // Host codepath

      //ComputeNeigh
      {
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPComputeNeighCPU>(chunk_size,team_size);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeNeighCPU> policy_neigh(chunk_size,team_size,vector_length);
        Kokkos::parallel_for("ComputeNeighCPU",policy_neigh,*this);
      }

      //PreUi
      {
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPPreUiCPU>(chunk_size,team_size);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUiCPU> policy_preui_cpu((chunk_size+team_size-1)/team_size,team_size,vector_length);
        Kokkos::parallel_for("PreUiCPU",policy_preui_cpu,*this);
      }

      // ComputeUi
      {
        int team_size = team_size_default;
        // Fused calculation of ulist and accumulation into ulisttot using atomics
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiCPU> policy_ui_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        Kokkos::parallel_for("ComputeUiCPU",policy_ui_cpu,*this);
      }

      {
        // Expand ulisttot -> ulisttot_full
        // Zero out ylist
        typename Kokkos::MDRangePolicy<DeviceType, Kokkos::IndexType<int>, Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>, TagPairSNAPTransformUiCPU> policy_transform_ui_cpu({0,0},{twojmax+1,chunk_size});
        Kokkos::parallel_for("TransformUiCPU",policy_transform_ui_cpu,*this);
      }

      //Compute bispectrum
      if (quadraticflag || eflag) {
        //ComputeZi
        int idxz_max = snaKK.idxz_max;
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeZiCPU> policy_zi_cpu(0,chunk_size*idxz_max);
        Kokkos::parallel_for("ComputeZiCPU",policy_zi_cpu,*this);

        //ComputeBi
        int team_size = team_size_default;
        check_team_size_for<TagPairSNAPComputeBiCPU>(chunk_size,team_size);
        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeBiCPU> policy_bi_cpu(chunk_size,team_size,vector_length);
        Kokkos::parallel_for("ComputeBiCPU",policy_bi_cpu,*this);
      }

      //ComputeYi
      {
        //Compute beta = dE_i/dB_i for all i in list
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPBetaCPU> policy_beta(0,chunk_size);
        Kokkos::parallel_for("ComputeBetaCPU",policy_beta,*this);

        //ComputeYi
        int idxz_max = snaKK.idxz_max;
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeYiCPU> policy_yi_cpu(0,chunk_size*idxz_max);
        Kokkos::parallel_for("ComputeYiCPU",policy_yi_cpu,*this);
      } // host flag

      //ComputeDuidrj and Deidrj
      {
        int team_size = team_size_default;

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDuidrjCPU> policy_duidrj_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);
        Kokkos::parallel_for("ComputeDuidrjCPU",policy_duidrj_cpu,*this);

        typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDeidrjCPU> policy_deidrj_cpu(((chunk_size+team_size-1)/team_size)*max_neighs,team_size,vector_length);

        Kokkos::parallel_for("ComputeDeidrjCPU",policy_deidrj_cpu,*this);
      }

    } else { // GPU

#ifdef LMP_KOKKOS_GPU

      // Pre-compute ceil(chunk_size / vector_length) for code cleanliness
      const int chunk_size_div = (chunk_size + vector_length - 1) / vector_length;

      //ComputeNeigh
      {
        // team_size_compute_neigh is defined in `pair_snap_kokkos.h`
        int scratch_size = scratch_size_helper<int>(team_size_compute_neigh * max_neighs);

        SnapAoSoATeamPolicy<DeviceType, team_size_compute_neigh, TagPairSNAPComputeNeigh> policy_neigh(chunk_size,team_size_compute_neigh,vector_length);
        policy_neigh = policy_neigh.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
        Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
      }

      //ComputeCayleyKlein
      {
        // tile_size_compute_ck is defined in `pair_snap_kokkos.h`
        Snap3DRangePolicy<DeviceType, tile_size_compute_ck, TagPairSNAPComputeCayleyKlein>
            policy_compute_ck({0,0,0},{vector_length,max_neighs,chunk_size_div},{vector_length,tile_size_compute_ck,1});
        Kokkos::parallel_for("ComputeCayleyKlein",policy_compute_ck,*this);
      }

      //PreUi
      {
        // tile_size_pre_ui is defined in `pair_snap_kokkos.h`
        Snap3DRangePolicy<DeviceType, tile_size_pre_ui, TagPairSNAPPreUi>
            policy_preui({0,0,0},{vector_length,twojmax+1,chunk_size_div},{vector_length,tile_size_pre_ui,1});
        Kokkos::parallel_for("PreUi",policy_preui,*this);
      }

      // ComputeUi w/vector parallelism, shared memory, direct atomicAdd into ulisttot
      {
        // team_size_compute_ui is defined in `pair_snap_kokkos.h`
        // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
        const int tile_size = vector_length * (twojmax + 1);
        const int scratch_size = scratch_size_helper<complex>(team_size_compute_ui * tile_size);

        if (chunk_size < parallel_thresh)
        {
          // Version with parallelism over j_bend

          // total number of teams needed: (natoms / 32) * (max_neighs) * ("bend" locations)
          const int n_teams = chunk_size_div * max_neighs * (twojmax + 1);
          const int n_teams_div = (n_teams + team_size_compute_ui - 1) / team_size_compute_ui;

          SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagPairSNAPComputeUiSmall> policy_ui(n_teams_div, team_size_compute_ui, vector_length);
          policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeUiSmall",policy_ui,*this);
        } else {
          // Version w/out parallelism over j_bend

          // total number of teams needed: (natoms / 32) * (max_neighs)
          const int n_teams = chunk_size_div * max_neighs;
          const int n_teams_div = (n_teams + team_size_compute_ui - 1) / team_size_compute_ui;

          SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagPairSNAPComputeUiLarge> policy_ui(n_teams_div, team_size_compute_ui, vector_length);
          policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeUiLarge",policy_ui,*this);
        }
      }

      //TransformUi: un-"fold" ulisttot, zero ylist
      {
        // team_size_transform_ui is defined in `pair_snap_kokkos.h`
        Snap3DRangePolicy<DeviceType, tile_size_transform_ui, TagPairSNAPTransformUi>
            policy_transform_ui({0,0,0},{vector_length,snaKK.idxu_max,chunk_size_div},{vector_length,tile_size_transform_ui,1});
        Kokkos::parallel_for("TransformUi",policy_transform_ui,*this);
      }

      //Compute bispectrum in AoSoA data layout, transform Bi
      if (quadraticflag || eflag) {
        // team_size_[compute_zi, compute_bi, transform_bi] are defined in `pair_snap_kokkos.h`

        //ComputeZi
        const int idxz_max = snaKK.idxz_max;
        Snap3DRangePolicy<DeviceType, tile_size_compute_zi, TagPairSNAPComputeZi>
            policy_compute_zi({0,0,0},{vector_length,idxz_max,chunk_size_div},{vector_length,tile_size_compute_zi,1});
        Kokkos::parallel_for("ComputeZi",policy_compute_zi,*this);

        //ComputeBi
        const int idxb_max = snaKK.idxb_max;
        Snap3DRangePolicy<DeviceType, tile_size_compute_bi, TagPairSNAPComputeBi>
            policy_compute_bi({0,0,0},{vector_length,idxb_max,chunk_size_div},{vector_length,tile_size_compute_bi,1});
        Kokkos::parallel_for("ComputeBi",policy_compute_bi,*this);

        //Transform data layout of blist out of AoSoA
        //We need this because `blist` gets used in ComputeForce which doesn't
        //take advantage of AoSoA, which at best would only be beneficial on the margins
        Snap3DRangePolicy<DeviceType, tile_size_transform_bi, TagPairSNAPTransformBi>
            policy_transform_bi({0,0,0},{vector_length,idxb_max,chunk_size_div},{vector_length,tile_size_transform_bi,1});
        Kokkos::parallel_for("TransformBi",policy_transform_bi,*this);
      }

      //Note zeroing `ylist` is fused into `TransformUi`.
      {
        //Compute beta = dE_i/dB_i for all i in list
        typename Kokkos::RangePolicy<DeviceType,TagPairSNAPBeta> policy_beta(0,chunk_size);
        Kokkos::parallel_for("ComputeBeta",policy_beta,*this);
        const int idxz_max = snaKK.idxz_max;
        if (quadraticflag || eflag) {
          Snap3DRangePolicy<DeviceType, tile_size_compute_yi, TagPairSNAPComputeYiWithZlist>
              policy_compute_yi({0,0,0},{vector_length,idxz_max,chunk_size_div},{vector_length,tile_size_compute_yi,1});
          Kokkos::parallel_for("ComputeYiWithZlist",policy_compute_yi,*this);
        } else {
          Snap3DRangePolicy<DeviceType, tile_size_compute_yi, TagPairSNAPComputeYi>
              policy_compute_yi({0,0,0},{vector_length,idxz_max,chunk_size_div},{vector_length,tile_size_compute_yi,1});
          Kokkos::parallel_for("ComputeYi",policy_compute_yi,*this);
        }
      }

      // Fused ComputeDuidrj, ComputeDeidrj
      {
        // team_size_compute_fused_deidrj is defined in `pair_snap_kokkos.h`

        // scratch size: 32 atoms * (twojmax+1) cached values * 2 for u, du, no double buffer
        const int tile_size = vector_length * (twojmax + 1);
        const int scratch_size = scratch_size_helper<complex>(2 * team_size_compute_fused_deidrj * tile_size);

        if (chunk_size < parallel_thresh)
        {
          // Version with parallelism over j_bend

          // total number of teams needed: (natoms / 32) * (max_neighs) * ("bend" locations)
          const int n_teams = chunk_size_div * max_neighs * (twojmax + 1);
          const int n_teams_div = (n_teams + team_size_compute_fused_deidrj - 1) / team_size_compute_fused_deidrj;

          // x direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjSmall<0> > policy_fused_deidrj_x(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_x = policy_fused_deidrj_x.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjSmall<0>",policy_fused_deidrj_x,*this);

          // y direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjSmall<1> > policy_fused_deidrj_y(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_y = policy_fused_deidrj_y.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjSmall<1>",policy_fused_deidrj_y,*this);

          // z direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjSmall<2> > policy_fused_deidrj_z(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_z = policy_fused_deidrj_z.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjSmall<2>",policy_fused_deidrj_z,*this);
        } else {
          // Version w/out parallelism over j_bend

          // total number of teams needed: (natoms / 32) * (max_neighs)
          const int n_teams = chunk_size_div * max_neighs;
          const int n_teams_div = (n_teams + team_size_compute_fused_deidrj - 1) / team_size_compute_fused_deidrj;

          // x direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjLarge<0> > policy_fused_deidrj_x(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_x = policy_fused_deidrj_x.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjLarge<0>",policy_fused_deidrj_x,*this);

          // y direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjLarge<1> > policy_fused_deidrj_y(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_y = policy_fused_deidrj_y.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjLarge<1>",policy_fused_deidrj_y,*this);

          // z direction
          SnapAoSoATeamPolicy<DeviceType, team_size_compute_fused_deidrj, TagPairSNAPComputeFusedDeidrjLarge<2> > policy_fused_deidrj_z(n_teams_div,team_size_compute_fused_deidrj,vector_length);
          policy_fused_deidrj_z = policy_fused_deidrj_z.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for("ComputeFusedDeidrjLarge<2>",policy_fused_deidrj_z,*this);

        }
      }

#endif // LMP_KOKKOS_GPU

    }

    //ComputeForce
    {
      if (evflag) {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeForce<HALF,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeForce<HALFTHREAD,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
        }
      } else {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeForce<HALF,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for(policy_force, *this);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairSNAPComputeForce<HALFTHREAD,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for(policy_force, *this);
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


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::allocate()
{
  PairSNAP::allocate();

  int n = atom->ntypes;
  MemKK::realloc_kokkos(d_map,"PairSNAPKokkos::map",n+1);

  MemKK::realloc_kokkos(k_cutsq,"PairSNAPKokkos::cutsq",n+1,n+1);
  rnd_cutsq = k_cutsq.template view<DeviceType>();
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
double PairSNAPKokkos<DeviceType, real_type, vector_length>::init_one(int i, int j)
{
  double cutone = PairSNAP::init_one(i,j);
  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::coeff(int narg, char **arg)
{
  PairSNAP::coeff(narg,arg);

  // Set up element lists

  MemKK::realloc_kokkos(d_radelem,"pair:radelem",nelements);
  MemKK::realloc_kokkos(d_wjelem,"pair:wjelem",nelements);
  MemKK::realloc_kokkos(d_coeffelem,"pair:coeffelem",nelements,ncoeffall);
  MemKK::realloc_kokkos(d_sinnerelem,"pair:sinnerelem",nelements);
  MemKK::realloc_kokkos(d_dinnerelem,"pair:dinnerelem",nelements);

  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_wjelem = Kokkos::create_mirror_view(d_wjelem);
  auto h_coeffelem = Kokkos::create_mirror_view(d_coeffelem);
  auto h_sinnerelem = Kokkos::create_mirror_view(d_sinnerelem);
  auto h_dinnerelem = Kokkos::create_mirror_view(d_dinnerelem);
  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int ielem = 0; ielem < nelements; ielem++) {
    h_radelem(ielem) = radelem[ielem];
    h_wjelem(ielem) = wjelem[ielem];
    h_sinnerelem(ielem) = sinnerelem[ielem];
    h_dinnerelem(ielem) = dinnerelem[ielem];
    for (int jcoeff = 0; jcoeff < ncoeffall; jcoeff++) {
      h_coeffelem(ielem,jcoeff) = coeffelem[ielem][jcoeff];
    }
  }

  for (int i = 1; i <= atom->ntypes; i++) {
    h_map(i) = map[i];
  }

  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_wjelem,h_wjelem);
  Kokkos::deep_copy(d_coeffelem,h_coeffelem);
  Kokkos::deep_copy(d_sinnerelem,h_sinnerelem);
  Kokkos::deep_copy(d_dinnerelem,h_dinnerelem);
  Kokkos::deep_copy(d_map,h_map);

  snaKK = SNAKokkos<DeviceType, real_type, vector_length>(rfac0,twojmax,
    rmin0,switchflag,bzeroflag,chemflag,bnormflag,wselfallflag,nelements,switchinnerflag);
  snaKK.grow_rij(0,0);
  snaKK.init();
}

/* ----------------------------------------------------------------------
   Begin routines that are unique to the GPU codepath. These take advantage
   of AoSoA data layouts and scratch memory for recursive polynomials
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPBeta,const int& ii) const {

  if (ii >= chunk_size) return;

  const int iatom_mod = ii % vector_length;
  const int iatom_div = ii / vector_length;

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const int ielem = d_map[itype];
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  auto d_coeffi = Kokkos::subview(d_coeffelem, ielem, Kokkos::ALL);

  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
    d_beta_pack(iatom_mod,icoeff,iatom_div) = d_coeffi[icoeff+1];
  }

  if (quadraticflag) {
    const auto idxb_max = my_sna.idxb_max;
    int k = ncoeff+1;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      const auto idxb = icoeff % idxb_max;
      const auto idx_chem = icoeff / idxb_max;
      real_type bveci = my_sna.blist(ii, idx_chem, idxb);
      d_beta_pack(iatom_mod,icoeff,iatom_div) += d_coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
        const auto jdxb = jcoeff % idxb_max;
        const auto jdx_chem = jcoeff / idxb_max;
        real_type bvecj = my_sna.blist(ii, jdx_chem, jdxb);
        d_beta_pack(iatom_mod,icoeff,iatom_div) += d_coeffi[k]*bvecj;
        d_beta_pack(iatom_mod,jcoeff,iatom_div) += d_coeffi[k]*bveci;
        k++;
      }
    }
  }
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeNeigh>::member_type& team) const {

  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // extract atom number
  int ii = team.team_rank() + team.league_rank() * team.team_size();
  if (ii >= chunk_size) return;

  // get a pointer to scratch memory
  // This is used to cache whether or not an atom is within the cutoff.
  // If it is, type_cache is assigned to the atom type.
  // If it's not, it's assigned to -1.
  const int tile_size = max_neighs; // number of elements per thread
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * tile_size; // offset into pointer for entire team
  int* type_cache = (int*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(int), 0) + scratch_shift;

  // Load various info about myself up front
  const int i = d_ilist[ii + chunk_offset];
  const F_FLOAT xtmp = x(i,0);
  const F_FLOAT ytmp = x(i,1);
  const F_FLOAT ztmp = x(i,2);
  const int itype = type[i];
  const int ielem = d_map[itype];
  const double radi = d_radelem[ielem];

  const int num_neighs = d_numneigh[i];

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // wj = weights for neighbors of I within cutoff
  // rcutij = cutoffs for neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  // Compute the number of neighbors, store rsq
  int ninside = 0;
  Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,num_neighs),
    [&] (const int jj, int& count) {
    T_INT j = d_neighbors(i,jj);
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;

    int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;

    if (rsq >= rnd_cutsq(itype,jtype)) {
      jtype = -1; // use -1 to signal it's outside the radius
    }

    type_cache[jj] = jtype;

    if (jtype >= 0)
     count++;
  }, ninside);

  d_ninside(ii) = ninside;

  Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_neighs),
    [&] (const int jj, int& offset, bool final) {

    const int jtype = type_cache[jj];

    if (jtype >= 0) {
      if (final) {
        T_INT j = d_neighbors(i,jj);
        const F_FLOAT dx = x(j,0) - xtmp;
        const F_FLOAT dy = x(j,1) - ytmp;
        const F_FLOAT dz = x(j,2) - ztmp;
        const int jelem = d_map[jtype];
        my_sna.rij(ii,offset,0) = static_cast<real_type>(dx);
        my_sna.rij(ii,offset,1) = static_cast<real_type>(dy);
        my_sna.rij(ii,offset,2) = static_cast<real_type>(dz);
        my_sna.wj(ii,offset) = static_cast<real_type>(d_wjelem[jelem]);
        my_sna.rcutij(ii,offset) = static_cast<real_type>((radi + d_radelem[jelem])*rcutfac);
        my_sna.inside(ii,offset) = j;
        if (switchinnerflag) {
          my_sna.sinnerij(ii,offset) = 0.5*(d_sinnerelem[ielem] + d_sinnerelem[jelem]);
          my_sna.dinnerij(ii,offset) = 0.5*(d_dinnerelem[ielem] + d_dinnerelem[jelem]);
        }
        if (chemflag)
          my_sna.element(ii,offset) = jelem;
        else
          my_sna.element(ii,offset) = 0;
      }
      offset++;
    }
  });
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeCayleyKlein,const int iatom_mod, const int jnbor, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int ii = iatom_mod + iatom_div * vector_length;
  if (ii >= chunk_size) return;

  const int ninside = d_ninside(ii);
  if (jnbor >= ninside) return;

  my_sna.compute_cayley_klein(iatom_mod,jnbor,iatom_div);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPPreUi, const int iatom_mod, const int j, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int ii = iatom_mod + iatom_div * vector_length;
  if (ii >= chunk_size) return;

  int itype = type(ii);
  int ielem = d_map[itype];

  my_sna.pre_ui(iatom_mod, j, ielem, iatom_div);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeUiSmall,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiSmall>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_ui;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / (max_neighs * (twojmax + 1)); // removed "const" to work around GCC 7 bug
  const int jj_jbend = flattened_idx - iatom_div * (max_neighs * (twojmax + 1));
  const int jbend = jj_jbend / max_neighs;
  int jj = jj_jbend - jbend * max_neighs; // removed "const" to work around GCC 7 bug

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    my_sna.compute_ui_small(team, iatom_mod, jbend, jj, iatom_div);
  });

}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeUiLarge,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiLarge>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_ui;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / max_neighs; // removed "const" to work around GCC 7 bug
  int jj = flattened_idx - iatom_div * max_neighs;

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    my_sna.compute_ui_large(team,iatom_mod, jj, iatom_div);
  });

}


template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPTransformUi,const int iatom_mod, const int idxu, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (idxu > my_sna.idxu_max) return;

  int elem_count = chemflag ? nelements : 1;

  for (int ielem = 0; ielem < elem_count; ielem++) {

    const FullHalfMapper mapper = my_sna.idxu_full_half[idxu];

    auto utot_re = my_sna.ulisttot_re_pack(iatom_mod, mapper.idxu_half, ielem, iatom_div);
    auto utot_im = my_sna.ulisttot_im_pack(iatom_mod, mapper.idxu_half, ielem, iatom_div);

    if (mapper.flip_sign == 1) {
      utot_im = -utot_im;
    } else if (mapper.flip_sign == -1) {
      utot_re = -utot_re;
    }

    my_sna.ulisttot_pack(iatom_mod, idxu, ielem, iatom_div) = { utot_re, utot_im };

    if (mapper.flip_sign == 0) {
      my_sna.ylist_pack_re(iatom_mod, mapper.idxu_half, ielem, iatom_div) = 0.;
      my_sna.ylist_pack_im(iatom_mod, mapper.idxu_half, ielem, iatom_div) = 0.;
    }
  }
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeYi,const int iatom_mod, const int jjz, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (jjz >= my_sna.idxz_max) return;

  my_sna.compute_yi(iatom_mod,jjz,iatom_div,d_beta_pack);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeYiWithZlist,const int iatom_mod, const int jjz, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (jjz >= my_sna.idxz_max) return;

  my_sna.compute_yi_with_zlist(iatom_mod,jjz,iatom_div,d_beta_pack);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeZi,const int iatom_mod, const int jjz, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (jjz >= my_sna.idxz_max) return;

  my_sna.compute_zi(iatom_mod,jjz,iatom_div);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeBi,const int iatom_mod, const int jjb, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (jjb >= my_sna.idxb_max) return;

  my_sna.compute_bi(iatom_mod,jjb,iatom_div);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPTransformBi,const int iatom_mod, const int idxb, const int iatom_div) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  if (idxb >= my_sna.idxb_max) return;

  const int ntriples = my_sna.ntriples;

  for (int itriple = 0; itriple < ntriples; itriple++) {

    const real_type blocal = my_sna.blist_pack(iatom_mod, idxb, itriple, iatom_div);

    my_sna.blist(iatom, itriple, idxb) = blocal;
  }

}

template<class DeviceType, typename real_type, int vector_length>
template<int dir>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeFusedDeidrjSmall<dir>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeFusedDeidrjSmall<dir> >::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_fused_deidrj;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / (max_neighs * (twojmax + 1)); // removed "const" to work around GCC 7 bug
  const int jj_jbend = flattened_idx - iatom_div * (max_neighs * (twojmax + 1));
  const int jbend = jj_jbend / max_neighs;
  int jj = jj_jbend - jbend * max_neighs; // removed "const" to work around GCC 7 bug

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    my_sna.template compute_fused_deidrj_small<dir>(team, iatom_mod, jbend, jj, iatom_div);

  });

}

template<class DeviceType, typename real_type, int vector_length>
template<int dir>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeFusedDeidrjLarge<dir>,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeFusedDeidrjLarge<dir> >::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_fused_deidrj;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / max_neighs; // removed "const" to work around GCC 7 bug
  int jj = flattened_idx - max_neighs * iatom_div;

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    my_sna.template compute_fused_deidrj_large<dir>(team, iatom_mod, jj, iatom_div);

  });
}

/* ----------------------------------------------------------------------
   Begin routines that are unique to the CPU codepath. These do not take
   advantage of AoSoA data layouts, but that could be a good point of
   future optimization and unification with the above kernels. It's unlikely
   that scratch memory optimizations will ever be useful for the CPU due to
   different arithmetic intensity requirements for the CPU vs GPU.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPBetaCPU,const int& ii) const {

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const int ielem = d_map[itype];
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  auto d_coeffi = Kokkos::subview(d_coeffelem, ielem, Kokkos::ALL);

  for (int icoeff = 0; icoeff < ncoeff; icoeff++)
    d_beta(icoeff,ii) = d_coeffi[icoeff+1];

  if (quadraticflag) {
    const auto idxb_max = my_sna.idxb_max;
    int k = ncoeff+1;
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      const auto idxb = icoeff % idxb_max;
      const auto idx_chem = icoeff / idxb_max;
      real_type bveci = my_sna.blist(ii,idx_chem,idxb);
      d_beta(icoeff,ii) += d_coeffi[k]*bveci;
      k++;
      for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
        const auto jdxb = jcoeff % idxb_max;
        const auto jdx_chem = jcoeff / idxb_max;
        real_type bvecj = my_sna.blist(ii,jdx_chem,jdxb);
        d_beta(icoeff,ii) += d_coeffi[k]*bvecj;
        d_beta(jcoeff,ii) += d_coeffi[k]*bveci;
        k++;
      }
    }
  }
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeNeighCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeNeighCPU>::member_type& team) const {


  int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;
  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);
  const int itype = type[i];
  const int ielem = d_map[itype];
  const double radi = d_radelem[ielem];

  const int num_neighs = d_numneigh[i];

  // rij[][3] = displacements between atom I and those neighbors
  // inside = indices of neighbors of I within cutoff
  // wj = weights for neighbors of I within cutoff
  // rcutij = cutoffs for neighbors of I within cutoff
  // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

  int ninside = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,num_neighs),
      [&] (const int jj, int& count) {
    Kokkos::single(Kokkos::PerThread(team), [&] () {
      T_INT j = d_neighbors(i,jj);
      const F_FLOAT dx = x(j,0) - xtmp;
      const F_FLOAT dy = x(j,1) - ytmp;
      const F_FLOAT dz = x(j,2) - ztmp;

      const int jtype = type(j);
      const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < rnd_cutsq(itype,jtype))
       count++;
    });
  },ninside);

  d_ninside(ii) = ninside;

  if (team.team_rank() == 0)
  Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team,num_neighs),
      [&] (const int jj, int& offset, bool final) {
  //for (int jj = 0; jj < num_neighs; jj++) {
    T_INT j = d_neighbors(i,jj);
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;

    const int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
    const int jelem = d_map[jtype];

    if (rsq < rnd_cutsq(itype,jtype)) {
      if (final) {
        my_sna.rij(ii,offset,0) = static_cast<real_type>(dx);
        my_sna.rij(ii,offset,1) = static_cast<real_type>(dy);
        my_sna.rij(ii,offset,2) = static_cast<real_type>(dz);
        my_sna.wj(ii,offset) = static_cast<real_type>(d_wjelem[jelem]);
        my_sna.rcutij(ii,offset) = static_cast<real_type>((radi + d_radelem[jelem])*rcutfac);
        my_sna.inside(ii,offset) = j;
        if (switchinnerflag) {
          my_sna.sinnerij(ii,offset) = 0.5*(d_sinnerelem[ielem] + d_sinnerelem[jelem]);
          my_sna.dinnerij(ii,offset) = 0.5*(d_dinnerelem[ielem] + d_dinnerelem[jelem]);
        }
        if (chemflag)
          my_sna.element(ii,offset) = jelem;
        else
          my_sna.element(ii,offset) = 0;
      }
      offset++;
    }
  });
}


template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPPreUiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPPreUiCPU>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // Extract the atom number
  const int ii = team.team_rank() + team.team_size() * team.league_rank();
  if (ii >= chunk_size) return;
  int itype = type(ii);
  int ielem = d_map[itype];

  my_sna.pre_ui_cpu(team,ii,ielem);
}



template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeUiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeUiCPU>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_ui_cpu(team,ii,jj);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPTransformUiCPU, const int j, const int iatom) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  if (iatom >= chunk_size) return;

  if (j > twojmax) return;

  int elem_count = chemflag ? nelements : 1;

  // De-symmetrize ulisttot
  for (int ielem = 0; ielem < elem_count; ielem++) {

    const int jju_half = my_sna.idxu_half_block(j);
    const int jju = my_sna.idxu_block(j);

    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        // Extract top half

        const int idxu_shift = mb * (j + 1) + ma;
        const int idxu_half = jju_half + idxu_shift;
        const int idxu = jju + idxu_shift;

        // Load ulist
        auto utot = my_sna.ulisttot(idxu_half, ielem, iatom);

        // Store
        my_sna.ulisttot_full(idxu, ielem, iatom) = utot;

        // Zero Yi
        my_sna.ylist(idxu_half, ielem, iatom) = {0., 0.};

        // Symmetric term
        const int sign_factor = (((ma+mb)%2==0)?1:-1);
        const int idxu_flip = jju + (j + 1 - mb) * (j + 1) - (ma + 1);

        if (sign_factor == 1) {
          utot.im = -utot.im;
        } else {
          utot.re = -utot.re;
        }

        my_sna.ulisttot_full(idxu_flip, ielem, iatom) = utot;
      }
    }
  }
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeYiCPU,const int& ii) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;
  my_sna.compute_yi_cpu(ii,d_beta);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeZiCPU,const int& ii) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;
  my_sna.compute_zi_cpu(ii);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeBiCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeBiCPU>::member_type& team) const {
  int ii = team.league_rank();
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;
  my_sna.compute_bi_cpu(team,ii);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeDuidrjCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDuidrjCPU>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_duidrj_cpu(team,ii,jj);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeDeidrjCPU,const typename Kokkos::TeamPolicy<DeviceType,TagPairSNAPComputeDeidrjCPU>::member_type& team) const {
  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() % ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ninside = d_ninside(ii);
  if (jj >= ninside) return;

  my_sna.compute_deidrj_cpu(team,ii,jj);
}

/* ----------------------------------------------------------------------
   Also used for both CPU and GPU codepaths. Could maybe benefit from a
   separate GPU/CPU codepath, but this kernel takes so little time it's
   likely not worth it.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>, const int& ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii + chunk_offset];

  SNAKokkos<DeviceType, real_type, vector_length> my_sna = snaKK;

  const int ninside = d_ninside(ii);

  for (int jj = 0; jj < ninside; jj++) {
    int j = my_sna.inside(ii,jj);

    F_FLOAT fij[3];
    fij[0] = my_sna.dedr(ii,jj,0);
    fij[1] = my_sna.dedr(ii,jj,1);
    fij[2] = my_sna.dedr(ii,jj,2);

    a_f(i,0) += fij[0];
    a_f(i,1) += fij[1];
    a_f(i,2) += fij[2];
    a_f(j,0) -= fij[0];
    a_f(j,1) -= fij[1];
    a_f(j,2) -= fij[2];
    // tally global and per-atom virial contribution
    if (EVFLAG) {
      if (vflag_either) {
        v_tally_xyz<NEIGHFLAG>(ev,i,j,
          fij[0],fij[1],fij[2],
          -my_sna.rij(ii,jj,0),-my_sna.rij(ii,jj,1),
          -my_sna.rij(ii,jj,2));
      }
    }

  }
  // tally energy contribution

  if (EVFLAG) {
    if (eflag_either) {

      const int itype = type(i);
      const int ielem = d_map[itype];
      auto d_coeffi = Kokkos::subview(d_coeffelem, ielem, Kokkos::ALL);

      // evdwl = energy of atom I, sum over coeffs_k * Bi_k

      auto evdwl = d_coeffi[0];

      // E = beta.B + 0.5*B^t.alpha.B

      const auto idxb_max = snaKK.idxb_max;

      // linear contributions

      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        const auto idxb = icoeff % idxb_max;
        const auto idx_chem = icoeff / idxb_max;
        evdwl += d_coeffi[icoeff+1]*my_sna.blist(ii,idx_chem,idxb);
      }

      // quadratic contributions
      if (quadraticflag) {
        int k = ncoeff+1;
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          const auto idxb = icoeff % idxb_max;
          const auto idx_chem = icoeff / idxb_max;
          real_type bveci = my_sna.blist(ii,idx_chem,idxb);
          evdwl += 0.5*d_coeffi[k++]*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            auto jdxb = jcoeff % idxb_max;
            auto jdx_chem = jcoeff / idxb_max;
            auto bvecj = my_sna.blist(ii,jdx_chem,jdxb);
            evdwl += d_coeffi[k++]*bveci*bvecj;
          }
        }
      }
      //ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
      if (eflag_global) ev.evdwl += evdwl;
      if (eflag_atom) d_eatom[i] += evdwl;
    }
  }
}

template<class DeviceType, typename real_type, int vector_length>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::operator() (TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSNAPComputeForce<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSNAPKokkos<DeviceType, real_type, vector_length>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
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

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
double PairSNAPKokkos<DeviceType, real_type, vector_length>::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += MemKK::memory_usage(d_beta);
  if (!host_flag)
    bytes += MemKK::memory_usage(d_beta_pack);
  bytes += MemKK::memory_usage(d_ninside);
  bytes += MemKK::memory_usage(d_map);
  bytes += MemKK::memory_usage(d_radelem);
  bytes += MemKK::memory_usage(d_wjelem);
  bytes += MemKK::memory_usage(d_coeffelem);
  bytes += MemKK::memory_usage(d_sinnerelem);
  bytes += MemKK::memory_usage(d_dinnerelem);
  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template<class TagStyle>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::check_team_size_for(int inum, int &team_size) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename real_type, int vector_length>
template<class TagStyle>
void PairSNAPKokkos<DeviceType, real_type, vector_length>::check_team_size_reduce(int inum, int &team_size) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename real_type, int vector_length>
template<typename scratch_type>
int PairSNAPKokkos<DeviceType, real_type, vector_length>::scratch_size_helper(int values_per_team) {
  typedef Kokkos::View<scratch_type*, Kokkos::DefaultExecutionSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}



/* ----------------------------------------------------------------------
   routines used by template reference classes
------------------------------------------------------------------------- */

template<class DeviceType>
PairSNAPKokkosDevice<DeviceType>::PairSNAPKokkosDevice(class LAMMPS *lmp)
   : PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>(lmp) { ; }

template<class DeviceType>
void PairSNAPKokkosDevice<DeviceType>::coeff(int narg, char **arg)
{
  Base::coeff(narg, arg);
}

template<class DeviceType>
void PairSNAPKokkosDevice<DeviceType>::init_style()
{
  Base::init_style();
}

template<class DeviceType>
double PairSNAPKokkosDevice<DeviceType>::init_one(int i, int j)
{
  return Base::init_one(i, j);
}

template<class DeviceType>
void PairSNAPKokkosDevice<DeviceType>::compute(int eflag_in, int vflag_in)
{
  Base::compute(eflag_in, vflag_in);
}

template<class DeviceType>
double PairSNAPKokkosDevice<DeviceType>::memory_usage()
{
  return Base::memory_usage();
}

#ifdef LMP_KOKKOS_GPU
template<class DeviceType>
PairSNAPKokkosHost<DeviceType>::PairSNAPKokkosHost(class LAMMPS *lmp)
   : PairSNAPKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>(lmp) { ; }

template<class DeviceType>
void PairSNAPKokkosHost<DeviceType>::coeff(int narg, char **arg)
{
  Base::coeff(narg, arg);
}

template<class DeviceType>
void PairSNAPKokkosHost<DeviceType>::init_style()
{
  Base::init_style();
}

template<class DeviceType>
double PairSNAPKokkosHost<DeviceType>::init_one(int i, int j)
{
  return Base::init_one(i, j);
}

template<class DeviceType>
void PairSNAPKokkosHost<DeviceType>::compute(int eflag_in, int vflag_in)
{
  Base::compute(eflag_in, vflag_in);
}

template<class DeviceType>
double PairSNAPKokkosHost<DeviceType>::memory_usage()
{
  return Base::memory_usage();
}
#endif

}
