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

#ifdef PAIR_CLASS

#else

// clang-format off
#ifndef LMP_PAIR_KOKKOS_H
#define LMP_PAIR_KOKKOS_H

#include "pair.h"               // IWYU pragma: export
#include "neighbor_kokkos.h"
#include "neigh_list_kokkos.h"
#include "math_special.h"
#include "update.h"
#include "Kokkos_Macros.hpp"
#include "Kokkos_ScatterView.hpp"

namespace LAMMPS_NS {

template<int Table>
struct CoulLongTable {
  enum {DoTable = Table};
};

// Tags for doing coulomb calculations or not
// They facilitate function overloading, since
// partial template specialization of member functions is not allowed
struct CoulTag {};
struct NoCoulTag {};

template<int FLAG>
struct DoCoul {
  typedef NoCoulTag type;
};

template<>
struct DoCoul<1> {
  typedef CoulTag type;
};


//Specialisation for Neighborlist types Half, HalfThread, Full
template <class PairStyle, int NEIGHFLAG, bool STACKPARAMS, int ZEROFLAG = 0, class Specialisation = void>
struct PairComputeFunctor  {
  typedef typename PairStyle::device_type device_type ;
  typedef ArrayTypes<device_type> AT;

  // Reduction type, contains evdwl, ecoul and virial[6]
  typedef EV_FLOAT value_type;

  // The copy of the pair style
  PairStyle c;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  int inum;

  using KKDeviceType = typename KKDevice<device_type>::value;
  using DUP = NeedDup_v<NEIGHFLAG,device_type>;

  // The force array is atomic for Half/Thread neighbor style
  //Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,
  //             typename KKDevice<device_type>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > f;
  KKScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,KKDeviceType,KKScatterSum,DUP> dup_f;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  //Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,
  //             typename KKDevice<device_type>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > eatom;
  KKScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,KKDeviceType,KKScatterSum,DUP> dup_eatom;

  //Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,
  //             typename KKDevice<device_type>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > vatom;
  KKScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,KKDeviceType,KKScatterSum,DUP> dup_vatom;

  NeighListKokkos<device_type> list;
  int use_cluster;   // 1 = use cluster-pair force kernel (GPU only)
  int nall;          // nlocal + nghost, for j-cluster bounds checking

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr),use_cluster(0),nall(0) {
    // allocate duplicated memory
    f = c.f;
    d_eatom = c.d_eatom;
    d_vatom = c.d_vatom;
    dup_f     = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(c.d_vatom);
    inum = list.inum;
  };

  // Set copymode = 1 so parent allocations aren't destructed by copies of the style
  ~PairComputeFunctor() {c.copymode = 1; list.copymode = 1;};

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }

  void contribute() {
    int need_dup = std::is_same_v<DUP,Kokkos::Experimental::ScatterDuplicated>;

    if (need_dup) {
      Kokkos::Experimental::contribute(c.f, dup_f);

      if (c.eflag_atom)
        Kokkos::Experimental::contribute(c.d_eatom, dup_eatom);

      if (c.vflag_atom)
        Kokkos::Experimental::contribute(c.d_vatom, dup_vatom);
    }
  }

  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const KK_FLOAT xtmp = c.x(i,0);
    const KK_FLOAT ytmp = c.x(i,1);
    const KK_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    KK_ACC_FLOAT fxtmp = 0;
    KK_ACC_FLOAT fytmp = 0;
    KK_ACC_FLOAT fztmp = 0;

    if (NEIGHFLAG == FULL && ZEROFLAG) {
      f(i,0) = 0;
      f(i,1) = 0;
      f(i,2) = 0;
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const KK_FLOAT delx = xtmp - c.x(j,0);
      const KK_FLOAT dely = ytmp - c.x(j,1);
      const KK_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        const KK_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

        fxtmp += static_cast<KK_ACC_FLOAT>(delx*fpair);
        fytmp += static_cast<KK_ACC_FLOAT>(dely*fpair);
        fztmp += static_cast<KK_ACC_FLOAT>(delz*fpair);

        if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || j < c.nlocal)) {
          a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
          a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
          a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);
        }

        if (EVFLAG) {
          KK_FLOAT evdwl = 0.0;
          if (c.eflag_either) {
            evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            const auto scale = (((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5));
            ev.evdwl += static_cast<KK_ACC_FLOAT>(scale *  evdwl);
          }

          if (c.vflag_either || c.eflag_atom) ev_tally(ev,i,j,evdwl,fpair,delx,dely,delz);
        }
      }

    }

    a_f(i,0) += static_cast<KK_ACC_FLOAT>(fxtmp);
    a_f(i,1) += static_cast<KK_ACC_FLOAT>(fytmp);
    a_f(i,2) += static_cast<KK_ACC_FLOAT>(fztmp);

    return ev;
  }

  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const KK_FLOAT xtmp = c.x(i,0);
    const KK_FLOAT ytmp = c.x(i,1);
    const KK_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);
    const KK_FLOAT qtmp = c.q(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    KK_ACC_FLOAT fxtmp = 0;
    KK_ACC_FLOAT fytmp = 0;
    KK_ACC_FLOAT fztmp = 0;

    if (NEIGHFLAG == FULL && ZEROFLAG) {
      f(i,0) = 0;
      f(i,1) = 0;
      f(i,2) = 0;
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
      const KK_FLOAT factor_coul = c.special_coul[sbmask(j)];
      j &= NEIGHMASK;
      const KK_FLOAT delx = xtmp - c.x(j,0);
      const KK_FLOAT dely = ytmp - c.x(j,1);
      const KK_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        KK_FLOAT fpair = KK_FLOAT();

        if (rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
          fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
        if (rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
          fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

        fxtmp += static_cast<KK_ACC_FLOAT>(delx*fpair);
        fytmp += static_cast<KK_ACC_FLOAT>(dely*fpair);
        fztmp += static_cast<KK_ACC_FLOAT>(delz*fpair);

        if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || j < c.nlocal)) {
          a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
          a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
          a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);
        }

        if (EVFLAG) {
          KK_FLOAT evdwl = 0.0;
          KK_FLOAT ecoul = 0.0;
          if (c.eflag_either) {
            if (rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype))) {
              evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              const auto scale = (((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5));
              ev.evdwl += static_cast<KK_ACC_FLOAT>(scale * evdwl);
            }
            if (rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype))) {
              ecoul = c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);
              const auto scale = (((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5));
              ev.ecoul += static_cast<KK_ACC_FLOAT>(scale * ecoul);
            }
          }

          if (c.vflag_either || c.eflag_atom) ev_tally(ev,i,j,evdwl+ecoul,fpair,delx,dely,delz);
        }
      }
    }

    a_f(i,0) += static_cast<KK_ACC_FLOAT>(fxtmp);
    a_f(i,1) += static_cast<KK_ACC_FLOAT>(fytmp);
    a_f(i,2) += static_cast<KK_ACC_FLOAT>(fztmp);

    return ev;
  }

  // TeamPolicy, newton off, and no energy/virial
  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel

  KOKKOS_FUNCTION
  void compute_item_team(typename Kokkos::TeamPolicy<device_type>::member_type team,
                         const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    const int atoms_per_team = team.team_size();
    const int firstatom = team.league_rank()*atoms_per_team;
    const int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const KK_FLOAT xtmp = c.x(i,0);
      const KK_FLOAT ytmp = c.x(i,1);
      const KK_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);

      if (NEIGHFLAG == FULL && ZEROFLAG) {
        Kokkos::single(Kokkos::PerThread(team), [&] (){
          f(i,0) = 0.0;
          f(i,1) = 0.0;
          f(i,2) = 0.0;
        });
      }

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      t_scalar3<KK_FLOAT> fsum;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, t_scalar3<KK_FLOAT>& ftmp) {

        int j = neighbors_i(jj);
        const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
        j &= NEIGHMASK;
        const KK_FLOAT delx = xtmp - c.x(j,0);
        const KK_FLOAT dely = ytmp - c.x(j,1);
        const KK_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          const KK_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

          const KK_FLOAT fx = delx*fpair;
          const KK_FLOAT fy = dely*fpair;
          const KK_FLOAT fz = delz*fpair;

          ftmp.x += fx;
          ftmp.y += fy;
          ftmp.z += fz;

          if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && j < c.nlocal) {
            a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
            a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
            a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
          }
        }

      },fsum);

      Kokkos::single(Kokkos::PerThread(team), [&] () {
        a_f(i,0) += static_cast<KK_ACC_FLOAT>(fsum.x);
        a_f(i,1) += static_cast<KK_ACC_FLOAT>(fsum.y);
        a_f(i,2) += static_cast<KK_ACC_FLOAT>(fsum.z);
      });

    });
  }

  // TeamPolicy, newton off, and no energy/virial
  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel

  KOKKOS_FUNCTION
  void compute_item_team(typename Kokkos::TeamPolicy<device_type>::member_type team,
                         const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    const int atoms_per_team = team.team_size();
    int firstatom = team.league_rank()*atoms_per_team;
    int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const KK_FLOAT xtmp = c.x(i,0);
      const KK_FLOAT ytmp = c.x(i,1);
      const KK_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);
      const KK_FLOAT qtmp = c.q(i);

      if (NEIGHFLAG == FULL && ZEROFLAG) {
        Kokkos::single(Kokkos::PerThread(team), [&] ()
        {
          f(i,0) = 0;
          f(i,1) = 0;
          f(i,2) = 0;
        });
      }

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      t_scalar3<KK_FLOAT> fsum;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, t_scalar3<KK_FLOAT>& ftmp) {
        int j = neighbors_i(jj);
        const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
        const KK_FLOAT factor_coul = c.special_coul[sbmask(j)];
        j &= NEIGHMASK;
        const KK_FLOAT delx = xtmp - c.x(j,0);
        const KK_FLOAT dely = ytmp - c.x(j,1);
        const KK_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          KK_FLOAT fpair = KK_FLOAT();

          if (rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
            fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          if (rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
            fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

          const KK_FLOAT fx = delx*fpair;
          const KK_FLOAT fy = dely*fpair;
          const KK_FLOAT fz = delz*fpair;

          ftmp.x += fx;
          ftmp.y += fy;
          ftmp.z += fz;

          if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && j < c.nlocal) {
            a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
            a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
            a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
          }
        }

      },fsum);

      Kokkos::single(Kokkos::PerThread(team), [&] () {
        a_f(i,0) += static_cast<KK_ACC_FLOAT>(fsum.x);
        a_f(i,1) += static_cast<KK_ACC_FLOAT>(fsum.y);
        a_f(i,2) += static_cast<KK_ACC_FLOAT>(fsum.z);
      });
    });
  }

  // TeamPolicy, newton off, and energy/virial
  // Loop over neighbors of one atom without coulomb interaction
  // This function is called in parallel

  KOKKOS_FUNCTION
  EV_FLOAT compute_item_team_ev(typename Kokkos::TeamPolicy<device_type>::member_type team,
                                const NeighListKokkos<device_type> &list, const NoCoulTag&) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;

    const int atoms_per_team = team.team_size();
    const int firstatom = team.league_rank()*atoms_per_team;
    const int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const KK_FLOAT xtmp = c.x(i,0);
      const KK_FLOAT ytmp = c.x(i,1);
      const KK_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);

      if (NEIGHFLAG == FULL && ZEROFLAG) {
        Kokkos::single(Kokkos::PerThread(team), [&] ()
        {
          f(i,0) = 0;
          f(i,1) = 0;
          f(i,2) = 0;
        });
      }

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      FEV_FLOAT fev;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, FEV_FLOAT& fev_tmp) {

        int j = neighbors_i(jj);
        const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
        j &= NEIGHMASK;
        const KK_FLOAT delx = xtmp - c.x(j,0);
        const KK_FLOAT dely = ytmp - c.x(j,1);
        const KK_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          const KK_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

          const KK_FLOAT fx = delx*fpair;
          const KK_FLOAT fy = dely*fpair;
          const KK_FLOAT fz = delz*fpair;

          fev_tmp.f[0] += static_cast<KK_ACC_FLOAT>(fx);
          fev_tmp.f[1] += static_cast<KK_ACC_FLOAT>(fy);
          fev_tmp.f[2] += static_cast<KK_ACC_FLOAT>(fz);

          const int I_CONTRIB = (NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD);
          const int J_CONTRIB = ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && j < c.nlocal);
          const KK_FLOAT factor = J_CONTRIB?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5);

          if (J_CONTRIB) {
            a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
            a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
            a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
          }

          KK_FLOAT evdwl = 0.0;
          if (c.eflag_either) {
            evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            fev_tmp.evdwl += static_cast<KK_ACC_FLOAT>(factor * evdwl);

            if (c.eflag_atom) {
              const KK_FLOAT epairhalf = static_cast<KK_FLOAT>(0.5) * evdwl;

              if (I_CONTRIB)
                a_eatom[i] += static_cast<KK_ACC_FLOAT>(epairhalf);

              if (J_CONTRIB)
                a_eatom[j] += static_cast<KK_ACC_FLOAT>(epairhalf);
            }
          }

          if (c.vflag_either) {
            const KK_FLOAT v_acc[6] = { delx*delx*fpair,
              dely*dely*fpair,
              delz*delz*fpair,
              delx*dely*fpair,
              delx*delz*fpair,
              dely*delz*fpair };

            const auto one_half = static_cast<KK_FLOAT>(0.5);

            for (int n = 0; n < 6; n++)
              fev_tmp.v[n] += static_cast<KK_ACC_FLOAT>(factor *v_acc[n]);

            if (c.vflag_atom) {
              if (I_CONTRIB) {
                for (int n = 0; n < 6; n++)
                  a_vatom(i, n) += static_cast<KK_ACC_FLOAT>(one_half * v_acc[n]);
              }
              if (J_CONTRIB) {
                for (int n = 0; n < 6; n++)
                  a_vatom(j, n) += static_cast<KK_ACC_FLOAT>(one_half * v_acc[n]);
              }
            }
          }
        }
      },fev);

      Kokkos::single(Kokkos::PerThread(team), [&] () {
        for (int n = 0; n < 3; n++)
          a_f(i,n) += static_cast<KK_ACC_FLOAT>(fev.f[n]);

        if (c.eflag_global)
          ev.evdwl += fev.evdwl;

        if (c.vflag_global) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += fev.v[n];
        }

        if (NEIGHFLAG == FULL) {

          if (c.eflag_atom)
            a_eatom(i) += fev.evdwl;

          if (c.vflag_atom) {
            for (int n = 0; n < 6; n++)
              a_vatom(i,n) += fev.v[n];
          }
        }
      });
    });
    return ev;
  }

  // TeamPolicy, newton off, and energy/virial
  // Loop over neighbors of one atom with coulomb interaction
  // This function is called in parallel

  KOKKOS_FUNCTION
  EV_FLOAT compute_item_team_ev(typename Kokkos::TeamPolicy<device_type>::member_type team,
                                const NeighListKokkos<device_type> &list, const CoulTag& ) const {

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;

    const int atoms_per_team = team.team_size();
    const int firstatom = team.league_rank()*atoms_per_team;
    const int lastatom = firstatom + atoms_per_team < inum ? firstatom + atoms_per_team : inum;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, firstatom, lastatom), [&] (const int &ii) {

      const int i = list.d_ilist[ii];
      const KK_FLOAT xtmp = c.x(i,0);
      const KK_FLOAT ytmp = c.x(i,1);
      const KK_FLOAT ztmp = c.x(i,2);
      const int itype = c.type(i);
      const KK_FLOAT qtmp = c.q(i);

      if (NEIGHFLAG == FULL && ZEROFLAG) {
        Kokkos::single(Kokkos::PerThread(team), [&] (){
          f(i,0) = 0.0;
          f(i,1) = 0.0;
          f(i,2) = 0.0;
        });
      }

      const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
      const int jnum = list.d_numneigh[i];

      FEV_FLOAT fev;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,jnum),
        [&] (const int jj, FEV_FLOAT& fev_tmp) {

        int j = neighbors_i(jj);
        const KK_FLOAT factor_lj = c.special_lj[sbmask(j)];
        const KK_FLOAT factor_coul = c.special_coul[sbmask(j)];
        j &= NEIGHMASK;
        const KK_FLOAT delx = xtmp - c.x(j,0);
        const KK_FLOAT dely = ytmp - c.x(j,1);
        const KK_FLOAT delz = ztmp - c.x(j,2);
        const int jtype = c.type(j);
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          KK_FLOAT fpair = KK_FLOAT();

          if (rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype)))
            fpair+=factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          if (rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype)))
            fpair+=c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);

          const KK_FLOAT fx = delx*fpair;
          const KK_FLOAT fy = dely*fpair;
          const KK_FLOAT fz = delz*fpair;

          fev_tmp.f[0] += static_cast<KK_ACC_FLOAT>(fx);
          fev_tmp.f[1] += static_cast<KK_ACC_FLOAT>(fy);
          fev_tmp.f[2] += static_cast<KK_ACC_FLOAT>(fz);

          const int I_CONTRIB = (NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD);
          const int J_CONTRIB = ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && j < c.nlocal);
          const KK_FLOAT factor = J_CONTRIB?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5);

          if (J_CONTRIB) {
            a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
            a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
            a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
          }

          KK_FLOAT evdwl = 0.0;
          KK_FLOAT ecoul = 0.0;
          if (c.eflag_either) {
            if (rsq < (STACKPARAMS?c.m_cut_ljsq[itype][jtype]:c.d_cut_ljsq(itype,jtype))) {
              evdwl = factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              fev_tmp.evdwl += static_cast<KK_ACC_FLOAT>(factor * evdwl);
            }
            if (rsq < (STACKPARAMS?c.m_cut_coulsq[itype][jtype]:c.d_cut_coulsq(itype,jtype))) {
              ecoul = c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,factor_coul,qtmp);
              fev_tmp.ecoul += static_cast<KK_ACC_FLOAT>(factor * ecoul);
            }


            if (c.eflag_atom) {
              const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * (evdwl + ecoul));

              if (I_CONTRIB)
                a_eatom[i] += epairhalf;

              if (J_CONTRIB)
                a_eatom[j] += epairhalf;
            }
          }

          if (c.vflag_either) {
            const KK_FLOAT v_acc[6] = { delx*delx*fpair,
              dely*dely*fpair,
              delz*delz*fpair,
              delx*dely*fpair,
              delx*delz*fpair,
              dely*delz*fpair };
            const auto one_half = static_cast<KK_FLOAT>(0.5);

            for (int n = 0; n < 6; n++)
              fev_tmp.v[n] += static_cast<KK_ACC_FLOAT>(factor * v_acc[n]);

            if (c.vflag_atom) {
              if (I_CONTRIB) {
                for (int n = 0; n < 6; n++)
                  a_vatom(i,n) += static_cast<KK_ACC_FLOAT>(one_half * v_acc[n]);
              }
              if (J_CONTRIB) {
                for (int n = 0; n < 6; n++)
                  a_vatom(j,n) += static_cast<KK_ACC_FLOAT>(one_half * v_acc[n]);
              }
            }
          }
        }
      },fev);

      Kokkos::single(Kokkos::PerThread(team), [&] () {
        for (int n = 0; n < 3; n++)
          a_f(i,n) += fev.f[n];

        if (c.eflag_global) {
          ev.evdwl += fev.evdwl;
          ev.ecoul += fev.ecoul;
        }

        if (c.vflag_global) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += fev.v[n];
        }

        if (NEIGHFLAG == FULL) {

          if (c.eflag_atom)
            a_eatom(i) += fev.evdwl + fev.ecoul;

          if (c.vflag_atom) {
            for (int n = 0; n < 6; n++)
              a_vatom(i,n) += fev.v[n];
          }
        }
      });
    });
    return ev;
  }

  // Cluster-pair kernel (GPU only, no special-pair scaling)
  // One warp (32 threads) per i-cluster of CI=8; loads CJ=4 j-coords into
  // scratch per tile and computes CI*CJ=32 pairs per tile (full warp utilization).
  // Supports FULL lists (newton off) and HALF/HALFTHREAD lists (newton on/off):
  // Phase D applies Newton 3rd law forces (and energy/virial) to j-atoms atomically.

  static constexpr int CI = 8;   // i-atoms per i-cluster (one warp = CI*CJ threads)
  static constexpr int CJ = 4;   // j-atoms per j-cluster tile

  // Scratch (non-EV path): 3*CI i-positions + 3*CJ j-positions + 3*CI force
  //   accumulators + 3*CI*CJ force tile + 2*CI ints (itype, iatom) + CJ ints (jtype)
  static constexpr int cluster_scratch_bytes =
      (6*CI + 3*CJ + 3*CI*CJ) * sizeof(KK_FLOAT) + (2*CI + CJ) * sizeof(int);
  // EV path additionally stores fpair and evdwl per tile for j-atom energy/virial
  static constexpr int cluster_scratch_bytes_ev =
      cluster_scratch_bytes + 2 * CI*CJ * sizeof(KK_FLOAT);

#if defined(LMP_KOKKOS_GPU)

  KOKKOS_FUNCTION
  void compute_item_cluster(
      typename Kokkos::TeamPolicy<device_type>::member_type team,
      const NeighListKokkos<device_type> &list, const NoCoulTag&) const
  {
    using ScratchF = Kokkos::View<KK_FLOAT*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;
    using ScratchI = Kokkos::View<int*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    const int ci = team.league_rank();
    const int start_ii = ci * CI;
    const int end_ii = (start_ii + CI < inum) ? start_ii + CI : inum;
    const int n_i = end_ii - start_ii;

    ScratchF s_xi(team.team_scratch(0), CI);
    ScratchF s_yi(team.team_scratch(0), CI);
    ScratchF s_zi(team.team_scratch(0), CI);
    ScratchF s_xj(team.team_scratch(0), CJ);
    ScratchF s_yj(team.team_scratch(0), CJ);
    ScratchF s_zj(team.team_scratch(0), CJ);
    ScratchF s_fxi(team.team_scratch(0), CI);
    ScratchF s_fyi(team.team_scratch(0), CI);
    ScratchF s_fzi(team.team_scratch(0), CI);
    // Force tile: one unique slot per thread — eliminates atomic_add contention
    ScratchF s_fx_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fy_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fz_tile(team.team_scratch(0), CI*CJ);
    ScratchI s_itype(team.team_scratch(0), CI);
    ScratchI s_jtype(team.team_scratch(0), CJ);
    ScratchI s_iatom(team.team_scratch(0), CI);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
      s_fxi(ki) = 0; s_fyi(ki) = 0; s_fzi(ki) = 0;
      if (ki < n_i) {
        const int i = list.d_ilist[start_ii + ki];
        s_iatom(ki) = i;
        s_xi(ki) = c.x(i,0); s_yi(ki) = c.x(i,1); s_zi(ki) = c.x(i,2);
        s_itype(ki) = c.type(i);
      }
    });
    team.team_barrier();

    const int cj_count = list.d_cluster_numneigh(ci);
    for (int cj_idx = 0; cj_idx < cj_count; cj_idx++) {
      const int cj = list.d_cluster_jlist(ci, cj_idx);
      const int first_j = cj * CJ;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CJ), [&](int kj) {
        const int j = first_j + kj;
        if (j < nall) {
          s_xj(kj) = c.x(j,0); s_yj(kj) = c.x(j,1); s_zj(kj) = c.x(j,2);
          s_jtype(kj) = c.type(j);
        } else {
          s_xj(kj) = static_cast<KK_FLOAT>(1e30);
          s_yj(kj) = static_cast<KK_FLOAT>(1e30);
          s_zj(kj) = static_cast<KK_FLOAT>(1e30);
          s_jtype(kj) = 1;
        }
      });
      team.team_barrier();

      // Phase B: CI*CJ = 32 active threads — full warp utilization.
      // Each thread writes to its unique tile slot; no atomic contention.
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ),
          [&](int pidx) {
        const int ki = pidx / CJ;
        const int kj = pidx % CJ;
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        if (ki < n_i) {
          const int i = s_iatom(ki);
          const int j = first_j + kj;
          if (i != j) {
            const KK_FLOAT delx = s_xi(ki) - s_xj(kj);
            const KK_FLOAT dely = s_yi(ki) - s_yj(kj);
            const KK_FLOAT delz = s_zi(ki) - s_zj(kj);
            const int itype = s_itype(ki);
            const int jtype = s_jtype(kj);
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < (STACKPARAMS ? c.m_cutsq[itype][jtype] : c.d_cutsq(itype,jtype))) {
              const KK_FLOAT fpair = c.template
                  compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              dfx = delx*fpair; dfy = dely*fpair; dfz = delz*fpair;
            }
          }
        }
        s_fx_tile(pidx) = dfx;
        s_fy_tile(pidx) = dfy;
        s_fz_tile(pidx) = dfz;
      });
      team.team_barrier();

      // Phase C: one thread per i-atom sums CJ tile slots — no atomics.
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        for (int kj = 0; kj < CJ; kj++) {
          dfx += s_fx_tile(ki*CJ + kj);
          dfy += s_fy_tile(ki*CJ + kj);
          dfz += s_fz_tile(ki*CJ + kj);
        }
        s_fxi(ki) += dfx;
        s_fyi(ki) += dfy;
        s_fzi(ki) += dfz;
      });
      team.team_barrier();

      // Phase D: Newton 3rd law — atomic force update to j-atoms (HALF/HALFTHREAD only).
      // For FULL lists each pair is computed from both sides; no Phase D needed.
      if (NEIGHFLAG != FULL) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ), [&](int pidx) {
          const int ki = pidx / CJ;
          const int kj = pidx % CJ;
          if (ki < n_i) {
            const int j = first_j + kj;
            if (j < nall && (c.newton_pair || j < c.nlocal)) {
              a_f(j, 0) -= static_cast<KK_ACC_FLOAT>(s_fx_tile(pidx));
              a_f(j, 1) -= static_cast<KK_ACC_FLOAT>(s_fy_tile(pidx));
              a_f(j, 2) -= static_cast<KK_ACC_FLOAT>(s_fz_tile(pidx));
            }
          }
        });
        team.team_barrier();
      }
    }

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n_i), [&](int ki) {
      a_f(s_iatom(ki),0) += static_cast<KK_ACC_FLOAT>(s_fxi(ki));
      a_f(s_iatom(ki),1) += static_cast<KK_ACC_FLOAT>(s_fyi(ki));
      a_f(s_iatom(ki),2) += static_cast<KK_ACC_FLOAT>(s_fzi(ki));
    });
  }

  KOKKOS_FUNCTION
  void compute_item_cluster(
      typename Kokkos::TeamPolicy<device_type>::member_type team,
      const NeighListKokkos<device_type> &list, const CoulTag&) const
  {
    using ScratchF = Kokkos::View<KK_FLOAT*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;
    using ScratchI = Kokkos::View<int*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    const int ci = team.league_rank();
    const int start_ii = ci * CI;
    const int end_ii = (start_ii + CI < inum) ? start_ii + CI : inum;
    const int n_i = end_ii - start_ii;

    ScratchF s_xi(team.team_scratch(0), CI);
    ScratchF s_yi(team.team_scratch(0), CI);
    ScratchF s_zi(team.team_scratch(0), CI);
    ScratchF s_xj(team.team_scratch(0), CJ);
    ScratchF s_yj(team.team_scratch(0), CJ);
    ScratchF s_zj(team.team_scratch(0), CJ);
    ScratchF s_fxi(team.team_scratch(0), CI);
    ScratchF s_fyi(team.team_scratch(0), CI);
    ScratchF s_fzi(team.team_scratch(0), CI);
    ScratchF s_fx_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fy_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fz_tile(team.team_scratch(0), CI*CJ);
    ScratchI s_itype(team.team_scratch(0), CI);
    ScratchI s_jtype(team.team_scratch(0), CJ);
    ScratchI s_iatom(team.team_scratch(0), CI);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
      s_fxi(ki) = 0; s_fyi(ki) = 0; s_fzi(ki) = 0;
      if (ki < n_i) {
        const int i = list.d_ilist[start_ii + ki];
        s_iatom(ki) = i;
        s_xi(ki) = c.x(i,0); s_yi(ki) = c.x(i,1); s_zi(ki) = c.x(i,2);
        s_itype(ki) = c.type(i);
      }
    });
    team.team_barrier();

    const int cj_count = list.d_cluster_numneigh(ci);
    for (int cj_idx = 0; cj_idx < cj_count; cj_idx++) {
      const int cj = list.d_cluster_jlist(ci, cj_idx);
      const int first_j = cj * CJ;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CJ), [&](int kj) {
        const int j = first_j + kj;
        if (j < nall) {
          s_xj(kj) = c.x(j,0); s_yj(kj) = c.x(j,1); s_zj(kj) = c.x(j,2);
          s_jtype(kj) = c.type(j);
        } else {
          s_xj(kj) = static_cast<KK_FLOAT>(1e30);
          s_yj(kj) = static_cast<KK_FLOAT>(1e30);
          s_zj(kj) = static_cast<KK_FLOAT>(1e30);
          s_jtype(kj) = 1;
        }
      });
      team.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ),
          [&](int pidx) {
        const int ki = pidx / CJ;
        const int kj = pidx % CJ;
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        if (ki < n_i) {
          const int i = s_iatom(ki);
          const int j = first_j + kj;
          if (i != j) {
            const KK_FLOAT delx = s_xi(ki) - s_xj(kj);
            const KK_FLOAT dely = s_yi(ki) - s_yj(kj);
            const KK_FLOAT delz = s_zi(ki) - s_zj(kj);
            const int itype = s_itype(ki);
            const int jtype = s_jtype(kj);
            const KK_FLOAT qtmp = c.q(i);
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < (STACKPARAMS ? c.m_cutsq[itype][jtype] : c.d_cutsq(itype,jtype))) {
              KK_FLOAT fpair = KK_FLOAT();
              if (rsq < (STACKPARAMS ? c.m_cut_ljsq[itype][jtype] : c.d_cut_ljsq(itype,jtype)))
                fpair += c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              if (rsq < (STACKPARAMS ? c.m_cut_coulsq[itype][jtype] : c.d_cut_coulsq(itype,jtype)))
                fpair += c.template compute_fcoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype,
                                                                               static_cast<KK_FLOAT>(1.0),qtmp);
              dfx = delx*fpair; dfy = dely*fpair; dfz = delz*fpair;
            }
          }
        }
        s_fx_tile(pidx) = dfx;
        s_fy_tile(pidx) = dfy;
        s_fz_tile(pidx) = dfz;
      });
      team.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        for (int kj = 0; kj < CJ; kj++) {
          dfx += s_fx_tile(ki*CJ + kj);
          dfy += s_fy_tile(ki*CJ + kj);
          dfz += s_fz_tile(ki*CJ + kj);
        }
        s_fxi(ki) += dfx;
        s_fyi(ki) += dfy;
        s_fzi(ki) += dfz;
      });
      team.team_barrier();

      if (NEIGHFLAG != FULL) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ), [&](int pidx) {
          const int ki = pidx / CJ;
          const int kj = pidx % CJ;
          if (ki < n_i) {
            const int j = first_j + kj;
            if (j < nall && (c.newton_pair || j < c.nlocal)) {
              a_f(j, 0) -= static_cast<KK_ACC_FLOAT>(s_fx_tile(pidx));
              a_f(j, 1) -= static_cast<KK_ACC_FLOAT>(s_fy_tile(pidx));
              a_f(j, 2) -= static_cast<KK_ACC_FLOAT>(s_fz_tile(pidx));
            }
          }
        });
        team.team_barrier();
      }
    }

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n_i), [&](int ki) {
      a_f(s_iatom(ki),0) += static_cast<KK_ACC_FLOAT>(s_fxi(ki));
      a_f(s_iatom(ki),1) += static_cast<KK_ACC_FLOAT>(s_fyi(ki));
      a_f(s_iatom(ki),2) += static_cast<KK_ACC_FLOAT>(s_fzi(ki));
    });
  }

  KOKKOS_FUNCTION
  EV_FLOAT compute_item_cluster_ev(
      typename Kokkos::TeamPolicy<device_type>::member_type team,
      const NeighListKokkos<device_type> &list, const NoCoulTag&) const
  {
    using ScratchF = Kokkos::View<KK_FLOAT*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;
    using ScratchI = Kokkos::View<int*,
        typename device_type::scratch_memory_space,
        Kokkos::MemoryUnmanaged>;

    auto a_f = dup_f.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    EV_FLOAT ev;

    const int ci = team.league_rank();
    const int start_ii = ci * CI;
    const int end_ii = (start_ii + CI < inum) ? start_ii + CI : inum;
    const int n_i = end_ii - start_ii;

    ScratchF s_xi(team.team_scratch(0), CI);
    ScratchF s_yi(team.team_scratch(0), CI);
    ScratchF s_zi(team.team_scratch(0), CI);
    ScratchF s_xj(team.team_scratch(0), CJ);
    ScratchF s_yj(team.team_scratch(0), CJ);
    ScratchF s_zj(team.team_scratch(0), CJ);
    ScratchF s_fxi(team.team_scratch(0), CI);
    ScratchF s_fyi(team.team_scratch(0), CI);
    ScratchF s_fzi(team.team_scratch(0), CI);
    ScratchF s_fx_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fy_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_fz_tile(team.team_scratch(0), CI*CJ);
    ScratchI s_itype(team.team_scratch(0), CI);
    ScratchI s_jtype(team.team_scratch(0), CJ);
    ScratchI s_iatom(team.team_scratch(0), CI);
    // Per-tile fpair and evdwl for j-atom energy/virial in Phase D (HALF lists only;
    // allocated unconditionally so scratch size is always cluster_scratch_bytes_ev)
    ScratchF s_fpair_tile(team.team_scratch(0), CI*CJ);
    ScratchF s_evdwl_tile(team.team_scratch(0), CI*CJ);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
      s_fxi(ki) = 0; s_fyi(ki) = 0; s_fzi(ki) = 0;
      if (ki < n_i) {
        const int i = list.d_ilist[start_ii + ki];
        s_iatom(ki) = i;
        s_xi(ki) = c.x(i,0); s_yi(ki) = c.x(i,1); s_zi(ki) = c.x(i,2);
        s_itype(ki) = c.type(i);
      }
    });
    team.team_barrier();

    const int cj_count = list.d_cluster_numneigh(ci);
    for (int cj_idx = 0; cj_idx < cj_count; cj_idx++) {
      const int cj = list.d_cluster_jlist(ci, cj_idx);
      const int first_j = cj * CJ;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CJ), [&](int kj) {
        const int j = first_j + kj;
        if (j < nall) {
          s_xj(kj) = c.x(j,0); s_yj(kj) = c.x(j,1); s_zj(kj) = c.x(j,2);
          s_jtype(kj) = c.type(j);
        } else {
          s_xj(kj) = static_cast<KK_FLOAT>(1e30);
          s_yj(kj) = static_cast<KK_FLOAT>(1e30);
          s_zj(kj) = static_cast<KK_FLOAT>(1e30);
          s_jtype(kj) = 1;
        }
      });
      team.team_barrier();

      // Phase B: forces, energy, virial per (ki,kj) pair; write to tile (no atomics).
      // Energy factor: FULL list counts each pair twice (0.5); HALF counts once (1.0).
      // Global virial: FULL = 1x, HALF+newton = 2x, HALF+no_newton = 1x (j's share
      //   added in Phase D via per-thread ev accumulator).
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ),
          [&](int pidx) {
        const int ki = pidx / CJ;
        const int kj = pidx % CJ;
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        KK_FLOAT fpair_v = 0, evdwl_v = 0;
        if (ki < n_i) {
          const int i = s_iatom(ki);
          const int j = first_j + kj;
          if (i != j) {
            const KK_FLOAT delx = s_xi(ki) - s_xj(kj);
            const KK_FLOAT dely = s_yi(ki) - s_yj(kj);
            const KK_FLOAT delz = s_zi(ki) - s_zj(kj);
            const int itype = s_itype(ki);
            const int jtype = s_jtype(kj);
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < (STACKPARAMS ? c.m_cutsq[itype][jtype] : c.d_cutsq(itype,jtype))) {
              const KK_FLOAT fpair = c.template
                  compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              fpair_v = fpair;
              dfx = delx*fpair; dfy = dely*fpair; dfz = delz*fpair;
              if (c.eflag_either || c.vflag_either) {
                const KK_FLOAT evdwl = c.template
                    compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
                evdwl_v = evdwl;
                const KK_FLOAT ef = (NEIGHFLAG == FULL) ?
                    static_cast<KK_FLOAT>(0.5) : static_cast<KK_FLOAT>(1.0);
                if (c.eflag_global)
                  ev.evdwl += static_cast<KK_ACC_FLOAT>(ef * evdwl);
                if (c.eflag_atom)
                  a_eatom[i] += static_cast<KK_ACC_FLOAT>(0.5*evdwl);
                if (c.vflag_global || c.vflag_atom) {
                  const auto half = static_cast<KK_FLOAT>(0.5);
                  const KK_FLOAT v[6] = { delx*delx*fpair*half, dely*dely*fpair*half,
                      delz*delz*fpair*half, delx*dely*fpair*half,
                      delx*delz*fpair*half, dely*delz*fpair*half };
                  if (c.vflag_global) {
                    if (NEIGHFLAG == FULL) {
                      for (int n = 0; n < 6; n++)
                        ev.v[n] += static_cast<KK_ACC_FLOAT>(v[n]);
                    } else {
                      // newton on: full pair virial (2x); newton off: i's share here,
                      // j's share added in Phase D when j is local
                      const KK_FLOAT vf = c.newton_pair ?
                          static_cast<KK_FLOAT>(2.0) : static_cast<KK_FLOAT>(1.0);
                      for (int n = 0; n < 6; n++)
                        ev.v[n] += static_cast<KK_ACC_FLOAT>(vf * v[n]);
                    }
                  }
                  if (c.vflag_atom)
                    for (int n = 0; n < 6; n++)
                      a_vatom(i,n) += static_cast<KK_ACC_FLOAT>(v[n]);
                }
              }
            }
          }
        }
        s_fx_tile(pidx) = dfx;
        s_fy_tile(pidx) = dfy;
        s_fz_tile(pidx) = dfz;
        s_fpair_tile(pidx) = fpair_v;
        s_evdwl_tile(pidx) = evdwl_v;
      });
      team.team_barrier();

      // Phase C: reduce tile into i-atom force accumulators — no atomics.
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI), [&](int ki) {
        KK_FLOAT dfx = 0, dfy = 0, dfz = 0;
        for (int kj = 0; kj < CJ; kj++) {
          dfx += s_fx_tile(ki*CJ + kj);
          dfy += s_fy_tile(ki*CJ + kj);
          dfz += s_fz_tile(ki*CJ + kj);
        }
        s_fxi(ki) += dfx;
        s_fyi(ki) += dfy;
        s_fzi(ki) += dfz;
      });
      team.team_barrier();

      // Phase D: Newton 3rd law — force, energy, and virial for j-atoms.
      if (NEIGHFLAG != FULL) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, CI*CJ), [&](int pidx) {
          const int ki = pidx / CJ;
          const int kj = pidx % CJ;
          if (ki < n_i) {
            const int j = first_j + kj;
            if (j < nall && (c.newton_pair || j < c.nlocal)) {
              a_f(j, 0) -= static_cast<KK_ACC_FLOAT>(s_fx_tile(pidx));
              a_f(j, 1) -= static_cast<KK_ACC_FLOAT>(s_fy_tile(pidx));
              a_f(j, 2) -= static_cast<KK_ACC_FLOAT>(s_fz_tile(pidx));
              if (c.eflag_atom)
                a_eatom[j] += static_cast<KK_ACC_FLOAT>(0.5 * s_evdwl_tile(pidx));
              if (c.vflag_atom || (!c.newton_pair && c.vflag_global)) {
                const KK_FLOAT fp = s_fpair_tile(pidx);
                const KK_FLOAT delx = s_xi(ki) - s_xj(kj);
                const KK_FLOAT dely = s_yi(ki) - s_yj(kj);
                const KK_FLOAT delz = s_zi(ki) - s_zj(kj);
                const auto half = static_cast<KK_FLOAT>(0.5);
                const KK_FLOAT v[6] = { delx*delx*fp*half, dely*dely*fp*half,
                    delz*delz*fp*half, delx*dely*fp*half,
                    delx*delz*fp*half, dely*delz*fp*half };
                if (c.vflag_atom)
                  for (int n = 0; n < 6; n++)
                    a_vatom(j,n) += static_cast<KK_ACC_FLOAT>(v[n]);
                // newton off: j's share of global virial (i < nlocal always in ilist)
                if (!c.newton_pair && c.vflag_global && j < c.nlocal)
                  for (int n = 0; n < 6; n++)
                    ev.v[n] += static_cast<KK_ACC_FLOAT>(v[n]);
              }
            }
          }
        });
        team.team_barrier();
      }
    }

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n_i), [&](int ki) {
      a_f(s_iatom(ki),0) += static_cast<KK_ACC_FLOAT>(s_fxi(ki));
      a_f(s_iatom(ki),1) += static_cast<KK_ACC_FLOAT>(s_fyi(ki));
      a_f(s_iatom(ki),2) += static_cast<KK_ACC_FLOAT>(s_fzi(ki));
    });

    return ev;
  }

  KOKKOS_FUNCTION
  EV_FLOAT compute_item_cluster_ev(
      typename Kokkos::TeamPolicy<device_type>::member_type team,
      const NeighListKokkos<device_type> &list, const CoulTag&) const
  {
    // Coul EV path: identical structure to NoCoul EV but adds Coulomb
    // Delegate to NoCoul EV for now (Coulomb EV can be added incrementally)
    return compute_item_cluster_ev(team, list, NoCoulTag{});
  }

#endif // LMP_KOKKOS_GPU

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const
  {
    auto a_eatom = dup_eatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();
    auto a_vatom = dup_vatom.template access<typename AtomicDup<NEIGHFLAG,device_type>::value>();

    const int EFLAG = c.eflag_either;
    const int NEWTON_PAIR = c.newton_pair;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * epair);
        if (NEWTON_PAIR || i < c.nlocal) a_eatom[i] += epairhalf;
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) a_eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const KK_FLOAT v0 = delx*delx*fpair;
      const KK_FLOAT v1 = dely*dely*fpair;
      const KK_FLOAT v2 = delz*delz*fpair;
      const KK_FLOAT v3 = delx*dely*fpair;
      const KK_FLOAT v4 = delx*delz*fpair;
      const KK_FLOAT v5 = dely*delz*fpair;
      const auto one_half = static_cast<KK_FLOAT>(0.5);

      const KK_ACC_FLOAT v_acc[6] = { static_cast<KK_ACC_FLOAT>(one_half*v0),
        static_cast<KK_ACC_FLOAT>(one_half*v1),
        static_cast<KK_ACC_FLOAT>(one_half*v2),
        static_cast<KK_ACC_FLOAT>(one_half*v3),
        static_cast<KK_ACC_FLOAT>(one_half*v4),
        static_cast<KK_ACC_FLOAT>(one_half*v5) };

      if (c.vflag_global) {
        if (NEIGHFLAG != FULL) {
          if (NEWTON_PAIR) {
            for (int n = 0; n < 6; n++)
              ev.v[n] += static_cast<KK_ACC_FLOAT>(2) * v_acc[n];
          } else {
            if (i < c.nlocal) {
              for (int n = 0; n < 6; n++)
                ev.v[n] += v_acc[n];
            }
            if (j < c.nlocal) {
              for (int n = 0; n < 6; n++)
                ev.v[n] += v_acc[n];
            }
          }
        } else {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_acc[n];
        }
      }

      if (c.vflag_atom) {
        if (NEWTON_PAIR || i < c.nlocal) {
          for (int n = 0; n < 6; n++)
            a_vatom(i,n) += v_acc[n];
        }
        if ((NEWTON_PAIR || j < c.nlocal) && NEIGHFLAG != FULL) {
          for (int n = 0; n < 6; n++)
            a_vatom(j,n) += v_acc[n];
        }
      }
    }
  }


// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (c.newton_pair) compute_item<0,1>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else compute_item<0,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    if (c.newton_pair)
      energy_virial += compute_item<1,1>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else
      energy_virial += compute_item<1,0>(i,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<device_type>::member_type& team) const {
#if defined(LMP_KOKKOS_GPU)
    if (use_cluster)
      compute_item_cluster(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else
#endif
      compute_item_team(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<device_type>::member_type& team, value_type &energy_virial) const {
#if defined(LMP_KOKKOS_GPU)
    if (use_cluster)
      energy_virial += compute_item_cluster_ev(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
    else
#endif
      energy_virial += compute_item_team_ev(team,list,typename DoCoul<PairStyle::COUL_FLAG>::type());
  }
};


// Filter out Neighflags which are not supported for PairStyle
// The enable_if clause will invalidate the last parameter of the function, so that
// a match is only achieved, if PairStyle supports the specific neighborlist variant.
// This uses the fact that failure to match template parameters is not an error.
// By having the enable_if with a ! and without it, exactly one of the functions
// pair_compute_neighlist will match - either the dummy version
// or the real one further below

template<class PairStyle, unsigned NEIGHFLAG, int ZEROFLAG = 0, class Specialisation = void>
EV_FLOAT pair_compute_neighlist (PairStyle* fpair, std::enable_if_t<!((NEIGHFLAG&PairStyle::EnabledNeighFlags) != 0), NeighListKokkos<typename PairStyle::device_type>*> list) {
  EV_FLOAT ev;
  (void) fpair;
  (void) list;
  printf("ERROR: calling pair_compute with invalid neighbor list style: requested %i  available %i \n",NEIGHFLAG,PairStyle::EnabledNeighFlags);
  return ev;
}

#if defined(LMP_KOKKOS_GPU)

// Build the cluster-pair neighbor list from the existing flat d_neighbors.
// Runs once per reneighbor as a post-pass.
// TeamPolicy: one warp (TEAM_SIZE threads) per i-cluster; dedup hash in team
// shared scratch (no DRAM spill).  All TEAM_SIZE threads cooperate to init,
// insert, and compact the hash.
//
// Shared scratch layout (int[HASH_SH + 2]):
//   [0 .. HASH_SH-1] : open-address hash table, sentinel = -1
//   [HASH_SH]        : unique j-cluster counter (populated in compact pass)
//   [HASH_SH+1]      : hash-overflow flag (1 if probe loop exhausted all slots)
//
// HASH_SH = 1024 -> 4 KB/team -> ~12 teams/SM on GP100 (48 KB shared).
// Handles up to ~800 unique j-clusters per i-cluster at a safe load factor.
// For larger cutoffs where unique j-clusters exceed this, increase HASH_SH
// and recompile (doubles shared mem usage per team, halving SM occupancy).

template<class DeviceType>
struct ClusterBuildFunctor {
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType> policy_type;
  typedef typename policy_type::member_type team_member;
  using ScratchI = Kokkos::View<int*, typename DeviceType::scratch_memory_space,
                                Kokkos::MemoryUnmanaged>;

  static constexpr int CI        = 8;    // i-atoms per i-cluster
  static constexpr int CJ        = 4;    // j-atoms per j-cluster
  static constexpr int TEAM_SIZE = 32;   // threads per team (one CUDA warp)
  static constexpr int HASH_SH   = 1024; // shared-mem hash slots

  typename AT::t_neighbors_2d_const d_neighbors;
  typename AT::t_int_1d_const       d_numneigh;
  typename AT::t_int_1d_const       d_ilist;
  typename AT::t_int_1d             d_cluster_numneigh;
  typename AT::t_int_2d             d_cluster_jlist;
  typename AT::t_int_1d             d_scratch;
  int inum;
  int max_jclusters;

  ClusterBuildFunctor(NeighListKokkos<DeviceType>* list) :
    d_neighbors(list->d_neighbors),
    d_numneigh(list->d_numneigh),
    d_ilist(list->d_ilist),
    d_cluster_numneigh(list->d_cluster_numneigh),
    d_cluster_jlist(list->d_cluster_jlist),
    d_scratch(list->d_cluster_scratch),
    inum(list->inum),
    max_jclusters(list->max_jclusters) {}

  static int scratch_size_needed() {
    return ScratchI::shmem_size(HASH_SH + 2);
  }

  KOKKOS_FUNCTION
  void operator()(const team_member& team) const {
    const int ci = team.league_rank();

    ScratchI scratch(team.team_scratch(0), HASH_SH + 2);

    // Init hash to sentinel, zero counter and overflow flag
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, HASH_SH), [&](int k) {
      scratch(k) = -1;
    });
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      scratch(HASH_SH)     = 0;
      scratch(HASH_SH + 1) = 0;
    });
    team.team_barrier();

    const int start_ii = ci * CI;
    const int end_ii   = (start_ii + CI < inum) ? start_ii + CI : inum;

    // Insert phase: distribute each i-atom's neighbor list across TEAM_SIZE
    // threads; each inserts its j-cluster into the shared hash via atomicCAS.
    // Single-warp team guarantees shared-mem coherence between iterations.
    for (int ii = start_ii; ii < end_ii; ii++) {
      if (scratch(HASH_SH + 1)) break; // hash saturated, no point continuing
      const int i    = d_ilist(ii);
      const int jnum = d_numneigh(i);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, jnum), [&](int jj) {
        const int cj = (d_neighbors(i, jj) & NEIGHMASK) / CJ;
        int slot = cj & (HASH_SH - 1);
        for (int probe = 0; probe < HASH_SH; probe++) {
          const int old = Kokkos::atomic_compare_exchange(&scratch(slot), -1, cj);
          if (old == -1 || old == cj) return; // inserted or already present
          slot = (slot + 1) & (HASH_SH - 1);
        }
        Kokkos::atomic_fetch_max(&scratch(HASH_SH + 1), 1); // hash full
      });
      team.team_barrier();
    }

    // Compact pass: scan hash in parallel, write unique j-clusters to output.
    // scratch(HASH_SH) starts at 0 (untouched during insert) and serves as
    // the output-slot counter; each thread atomically claims its write index.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, HASH_SH), [&](int k) {
      if (scratch(k) != -1) {
        const int idx = Kokkos::atomic_fetch_add(&scratch(HASH_SH), 1);
        if (idx < max_jclusters) d_cluster_jlist(ci, idx) = scratch(k);
      }
    });
    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      const int njc        = scratch(HASH_SH);
      const int hash_ovf   = scratch(HASH_SH + 1);
      const bool jlist_ovf = (njc > max_jclusters);
      d_cluster_numneigh(ci) = (hash_ovf || jlist_ovf) ? 0 : njc;
      if (hash_ovf)
        Kokkos::atomic_fetch_max(&d_scratch(0), 2);
      if (jlist_ovf) {
        Kokkos::atomic_fetch_max(&d_scratch(0), 1);
        Kokkos::atomic_fetch_max(&d_scratch(1), njc + 64);
      }
    });
  }
};

template<class DeviceType>
void build_cluster_list(NeighListKokkos<DeviceType>* list)
{
  if (list->inum == 0) return;
  const int CI            = ClusterBuildFunctor<DeviceType>::CI;
  const int TEAM_SIZE     = ClusterBuildFunctor<DeviceType>::TEAM_SIZE;
  const int num_iclusters = (list->inum + CI - 1) / CI;
  const int scratch_bytes = ClusterBuildFunctor<DeviceType>::scratch_size_needed();

  if (list->max_jclusters == 0) list->max_jclusters = 256;

  using PolicyType = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;

  for (;;) {
    list->grow_clusters(num_iclusters, list->max_jclusters);

    auto h = Kokkos::create_mirror_view(list->d_cluster_scratch);
    h(0) = 0; h(1) = list->max_jclusters; h(2) = 0;
    Kokkos::deep_copy(list->d_cluster_scratch, h);

    ClusterBuildFunctor<DeviceType> ff(list);
    PolicyType policy(num_iclusters, TEAM_SIZE);
    policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));
    Kokkos::parallel_for(policy, ff);
    Kokkos::fence();

    auto h2 = Kokkos::create_mirror_view(list->d_cluster_scratch);
    Kokkos::deep_copy(h2, list->d_cluster_scratch);
    if (h2(0) == 0) break;
    if (h2(0) >= 2)
      list->cluster_fatal(FLERR,
          "ClusterBuildFunctor: shared-mem hash saturated (" +
          std::to_string(ClusterBuildFunctor<DeviceType>::HASH_SH) +
          " slots). Increase HASH_SH in pair_kokkos.h and recompile, or reduce cutoff.");
    list->max_jclusters = h2(1);
  }
}

#endif // LMP_KOKKOS_GPU

template<class NeighStyle>
int GetMaxNeighs(NeighStyle* list)
{
  auto d_ilist = list->d_ilist;
  auto d_numneigh = list->d_numneigh;
  int inum = list->inum;

  int maxneigh = 0;
  Kokkos::parallel_reduce(inum, LAMMPS_LAMBDA(const int ii, int &maxneigh) {
    const int i = d_ilist[ii];
    const int num_neighs = d_numneigh[i];
    maxneigh = MAX(maxneigh,num_neighs);
  }, Kokkos::Max<int>(maxneigh));

  if (maxneigh < 0) maxneigh = 0;

  return maxneigh;
}

template<class DeviceType, class FunctorStyle>
void GetMaxTeamSize(FunctorStyle& functor, int inum,
                int &teamsize_max_for, int &teamsize_max_reduce)
{
  teamsize_max_for = Kokkos::TeamPolicy<DeviceType>(inum,Kokkos::AUTO).team_size_max(functor,Kokkos::ParallelForTag());
  teamsize_max_reduce = Kokkos::TeamPolicy<DeviceType>(inum,Kokkos::AUTO).team_size_max(functor,Kokkos::ParallelReduceTag());
}

// Submit ParallelFor for NEIGHFLAG=HALF,HALFTHREAD,FULL
template<class PairStyle, unsigned NEIGHFLAG, int ZEROFLAG = 0, class Specialisation = void>
EV_FLOAT pair_compute_neighlist (PairStyle* fpair, std::enable_if_t<(NEIGHFLAG&PairStyle::EnabledNeighFlags) != 0, NeighListKokkos<typename PairStyle::device_type>*> list) {
  EV_FLOAT ev;

  const int inum = list->inum;

  if (!fpair->lmp->kokkos->neigh_thread_set)
    if (fpair->lmp->kokkos->ngpus && inum <= 16000)
      if (NEIGHFLAG == FULL || !fpair->newton_pair)
        fpair->lmp->kokkos->neigh_thread = 1;

  if (fpair->lmp->kokkos->neigh_thread) {

    static int vectorsize = 0;
    static int atoms_per_team = 0;

#if defined(LMP_KOKKOS_GPU)
    static int teamsize_max_for = 0;
    static int teamsize_max_reduce = 0;
    static int lastcall = -1;
    if (!vectorsize || lastcall < fpair->lmp->neighbor->lastcall) {
      lastcall = fpair->lmp->update->ntimestep;
      vectorsize = GetMaxNeighs(list);
      if (vectorsize == 0) vectorsize = 1;
      vectorsize = static_cast<int>(MathSpecial::powint(2.0,(int(log2(double(vectorsize)) + 0.5)))); // round to nearest power of 2

  #if defined(KOKKOS_ENABLE_HIP)
      int max_vectorsize = 64;
  #else
      int max_vectorsize = 32;
  #endif

      if (fpair->lmp->kokkos->threads_per_atom_set)
        vectorsize = fpair->lmp->kokkos->threads_per_atom;

      vectorsize = MIN(vectorsize,max_vectorsize);

      if (fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
        PairComputeFunctor<PairStyle,NEIGHFLAG,false,ZEROFLAG,Specialisation > ff(fpair,list);
        GetMaxTeamSize<typename PairStyle::device_type>(ff, inum, teamsize_max_for, teamsize_max_reduce);
      } else {
        PairComputeFunctor<PairStyle,NEIGHFLAG,true,ZEROFLAG,Specialisation > ff(fpair,list);
        GetMaxTeamSize<typename PairStyle::device_type>(ff, inum, teamsize_max_for, teamsize_max_reduce);
      }

      // Build cluster list when neigh/cluster is enabled
      if (fpair->lmp->kokkos->neigh_cluster &&
          fpair->special_lj[1] == 1.0 && fpair->special_lj[2] == 1.0 &&
          fpair->special_lj[3] == 1.0)
        build_cluster_list<typename PairStyle::device_type>(list);
    }

    // Cluster-pair kernel: 32 threads per i-cluster, no vector lanes.
    // Supports FULL, HALF, and HALFTHREAD neighbor lists.
    const bool do_cluster = fpair->lmp->kokkos->neigh_cluster &&
        fpair->special_lj[1] == 1.0 && fpair->special_lj[2] == 1.0 &&
        fpair->special_lj[3] == 1.0;

    if (do_cluster) {
      using DeviceType = typename PairStyle::device_type;
      using PolicyType = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;

      using PCF = PairComputeFunctor<PairStyle,NEIGHFLAG,true,ZEROFLAG,Specialisation>;
      constexpr int cluster_ts = PCF::CI * PCF::CJ;  // CI*CJ = 32 threads per team
      const int num_iclusters = (inum + PCF::CI - 1) / PCF::CI;
      const int scratch_bytes = (fpair->eflag || fpair->vflag) ?
          PCF::cluster_scratch_bytes_ev : PCF::cluster_scratch_bytes;

      const int nall = fpair->atom->nlocal + fpair->atom->nghost;
      if (fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
        PairComputeFunctor<PairStyle,NEIGHFLAG,false,ZEROFLAG,Specialisation> ff(fpair,list);
        ff.use_cluster = 1; ff.nall = nall;
        PolicyType policy(num_iclusters, cluster_ts, 1);
        policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));
        if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
        else                              Kokkos::parallel_for(policy,ff);
        ff.contribute();
      } else {
        PairComputeFunctor<PairStyle,NEIGHFLAG,true,ZEROFLAG,Specialisation> ff(fpair,list);
        ff.use_cluster = 1; ff.nall = nall;
        PolicyType policy(num_iclusters, cluster_ts, 1);
        policy = policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));
        if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
        else                              Kokkos::parallel_for(policy,ff);
        ff.contribute();
      }
      return ev;
    }

    int teamsize_max = teamsize_max_for;
    if (fpair->eflag || fpair->vflag)
      teamsize_max = teamsize_max_reduce;

    if (fpair->lmp->kokkos->pair_team_size_set)
      teamsize_max = fpair->lmp->kokkos->pair_team_size;

    atoms_per_team = teamsize_max/vectorsize;
#else
    vectorsize = 1;
    atoms_per_team = 1;
#endif

    const int num_teams = inum / atoms_per_team + (inum % atoms_per_team ? 1 : 0);

    if (fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
      PairComputeFunctor<PairStyle,NEIGHFLAG,false,ZEROFLAG,Specialisation > ff(fpair,list);
      Kokkos::TeamPolicy<typename PairStyle::device_type,Kokkos::IndexType<int> > policy(num_teams,atoms_per_team,vectorsize);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
      else                              Kokkos::parallel_for(policy,ff);
      ff.contribute();
    } else {
      PairComputeFunctor<PairStyle,NEIGHFLAG,true,ZEROFLAG,Specialisation > ff(fpair,list);
      Kokkos::TeamPolicy<typename PairStyle::device_type,Kokkos::IndexType<int> > policy(num_teams,atoms_per_team,vectorsize);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(policy,ff,ev);
      else                              Kokkos::parallel_for(policy,ff);
      ff.contribute();
    }
  } else {
    if (fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
      PairComputeFunctor<PairStyle,NEIGHFLAG,false,ZEROFLAG,Specialisation > ff(fpair,list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(inum,ff,ev);
      else                              Kokkos::parallel_for(inum,ff);
      ff.contribute();
    } else {
      PairComputeFunctor<PairStyle,NEIGHFLAG,true,ZEROFLAG,Specialisation > ff(fpair,list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(inum,ff,ev);
      else                              Kokkos::parallel_for(inum,ff);
      ff.contribute();
    }
  }
  return ev;
}

template<class PairStyle, class Specialisation = void>
EV_FLOAT pair_compute (PairStyle* fpair, NeighListKokkos<typename PairStyle::device_type>* list) {
  EV_FLOAT ev;
  if (fpair->neighflag == FULL) {
    if (utils::strmatch(fpair->lmp->force->pair_style,"^hybrid")) {
      fpair->fuse_force_clear_flag = 0;
      ev = pair_compute_neighlist<PairStyle,FULL,0,Specialisation> (fpair,list);
    } else {
      fpair->fuse_force_clear_flag = 1;
      ev = pair_compute_neighlist<PairStyle,FULL,1,Specialisation> (fpair,list);
    }
  } else if (fpair->neighflag == HALFTHREAD) {
    ev = pair_compute_neighlist<PairStyle,HALFTHREAD,0,Specialisation> (fpair,list);
  } else if (fpair->neighflag == HALF) {
    ev = pair_compute_neighlist<PairStyle,HALF,0,Specialisation> (fpair,list);
  }
  return ev;
}

template<class DeviceType>
struct PairVirialFDotRCompute {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;
  typename AT::t_kkfloat_1d_3_lr_const_um x;
  typename AT::t_kkacc_1d_3_const_um f;
  const int offset;

  PairVirialFDotRCompute(  typename AT::t_kkfloat_1d_3_lr_const_um x_,
  typename AT::t_kkacc_1d_3_const_um f_,
  const int offset_):x(x_),f(f_),offset(offset_) {}

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, value_type &energy_virial) const {
    const int i = j + offset;
    energy_virial.v[0] += f(i,0)*static_cast<KK_ACC_FLOAT>(x(i,0));
    energy_virial.v[1] += f(i,1)*static_cast<KK_ACC_FLOAT>(x(i,1));
    energy_virial.v[2] += f(i,2)*static_cast<KK_ACC_FLOAT>(x(i,2));
    energy_virial.v[3] += f(i,1)*static_cast<KK_ACC_FLOAT>(x(i,0));
    energy_virial.v[4] += f(i,2)*static_cast<KK_ACC_FLOAT>(x(i,0));
    energy_virial.v[5] += f(i,2)*static_cast<KK_ACC_FLOAT>(x(i,1));
  }
};

template<class PairStyle>
void pair_virial_fdotr_compute(PairStyle* fpair) {
  EV_FLOAT virial;
  if (fpair->neighbor->includegroup == 0) {
    int nall = fpair->atom->nlocal + fpair->atom->nghost;
    Kokkos::parallel_reduce(nall,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,0),virial);
  } else {
    Kokkos::parallel_reduce(fpair->atom->nfirst,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,0),virial);
    EV_FLOAT virial_ghost;
    Kokkos::parallel_reduce(fpair->atom->nghost,PairVirialFDotRCompute<typename PairStyle::device_type>(fpair->x,fpair->f,fpair->atom->nlocal),virial_ghost);
    virial+=virial_ghost;
  }
  fpair->vflag_fdotr = 0;
  for (int n = 0; n < 6; n++)
    fpair->virial[n] = static_cast<double>(virial.v[n]);
}

}

#endif
#endif
