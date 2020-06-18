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
   Contributing author: Dan Ibanez (SNL)
------------------------------------------------------------------------- */

#include "pair_table_rx_kokkos.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "kokkos.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "fix.h"
#include "kokkos_few.h"
#include "kokkos.h"
#include "modify.h"
#include "utils.h"
#include <cassert>

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define OneFluidValue (-1)
#define isOneFluid(_site_) ( (_site_) == OneFluidValue )

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void getMixingWeights(
    typename ArrayTypes<Space>::t_float_2d_randomread dvector,
    int nspecies,
    int isite1, int isite2,
    bool fractionalWeighting,
    int id,
    typename GetFloatType<Space>::type &mixWtSite1old, typename GetFloatType<Space>::type &mixWtSite2old,
    typename GetFloatType<Space>::type &mixWtSite1, typename GetFloatType<Space>::type &mixWtSite2) {
  KK_FLOAT fractionOFAold, fractionOFA;
  KK_FLOAT fractionOld1, fraction1;
  KK_FLOAT fractionOld2, fraction2;
  KK_FLOAT nMoleculesOFAold, nMoleculesOFA;
  KK_FLOAT nMoleculesOld1, nMolecules1;
  KK_FLOAT nMoleculesOld2, nMolecules2;
  KK_FLOAT nTotal, nTotalOld;

  nTotal = 0.0;
  nTotalOld = 0.0;
  assert(id >= 0);
  assert(id < dvector.extent(1));
  for (int ispecies = 0; ispecies < nspecies; ++ispecies){
    assert(ispecies < dvector.extent(0));
    nTotal += dvector(ispecies,id);
    assert(ispecies+nspecies < dvector.extent(0));
    nTotalOld += dvector(ispecies+nspecies,id);
  }

  assert(isite1 >= 0);
  assert(isite1 < nspecies);
  assert(isite2 >= 0);
  assert(isite2 < nspecies);
  if (isOneFluid(isite1) == false){
    nMoleculesOld1 = dvector(isite1+nspecies,id);
    nMolecules1 = dvector(isite1,id);
    fractionOld1 = nMoleculesOld1/nTotalOld;
    fraction1 = nMolecules1/nTotal;
  }
  if (isOneFluid(isite2) == false){
    nMoleculesOld2 = dvector(isite2+nspecies,id);
    nMolecules2 = dvector(isite2,id);
    fractionOld2 = nMoleculesOld2/nTotalOld;
    fraction2 = nMolecules2/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)){
    nMoleculesOFAold  = 0.0;
    nMoleculesOFA  = 0.0;
    fractionOFAold  = 0.0;
    fractionOFA  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      if (isite1 == ispecies || isite2 == ispecies) continue;
      nMoleculesOFAold += dvector(ispecies+nspecies,id);
      nMoleculesOFA += dvector(ispecies,id);
      fractionOFAold += dvector(ispecies+nspecies,id)/nTotalOld;
      fractionOFA += dvector(ispecies,id)/nTotal;
    }
    if(isOneFluid(isite1)){
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if(isOneFluid(isite2)){
      nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld2 = fractionOFAold;
      fraction2 = fractionOFA;
    }
  }

  if(fractionalWeighting){
    mixWtSite1old = fractionOld1;
    mixWtSite1 = fraction1;
    mixWtSite2old = fractionOld2;
    mixWtSite2 = fraction2;
  } else {
    mixWtSite1old = nMoleculesOld1;
    mixWtSite1 = nMolecules1;
    mixWtSite2old = nMoleculesOld2;
    mixWtSite2 = nMolecules2;
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairTableRXKokkos<Space>::PairTableRXKokkos(LAMMPS *lmp) : PairTable(lmp)
{
  update_table = 0;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK |
                  DVECTOR_MASK | UCG_MASK | UCGNEW_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK;
  k_table = new TableDual();
  h_table = new TableHost();
  d_table = new TableDevice();
  fractionalWeighting = true;

  site1 = nullptr;
  site2 = nullptr;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairTableRXKokkos<Space>::~PairTableRXKokkos()
{
  if (copymode) return;

  delete [] site1;
  delete [] site2;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  if (allocated) {
    memoryKK->destroy_kokkos(k_table->k_cutsq, cutsq);
    memoryKK->destroy_kokkos(k_table->k_tabindex, tabindex);
  }

  delete k_table;
  k_table = nullptr;
  delete h_table;
  h_table = nullptr;
  delete d_table;
  d_table = nullptr;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  if(update_table)
    create_kokkos_tables();
  if(tabstyle == LOOKUP)
    compute_style<LOOKUP>(eflag_in,vflag_in);
  if(tabstyle == LINEAR)
    compute_style<LINEAR>(eflag_in,vflag_in);
  if(tabstyle == SPLINE)
    compute_style<SPLINE>(eflag_in,vflag_in);
  if(tabstyle == BITMAP)
    compute_style<BITMAP>(eflag_in,vflag_in);

  copymode = 0;
}

KOKKOS_INLINE_FUNCTION static int sbmask(const int& j)
{
  return j >> SBBITS & 3;
}

template <ExecutionSpace Space, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
static KK_FLOAT
compute_fpair(KK_FLOAT rsq,
              int itype, int jtype,
              typename PairTableRXKokkos<Space>::TableDeviceConst const& d_table_const
              ) {
  Pair::union_int_float_t rsq_lookup;
  KK_FLOAT fpair;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (TABSTYLE == PairTable::LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    fpair = d_table_const.f(tidx,itable);
  } else if (TABSTYLE == PairTable::LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  } else if (TABSTYLE == PairTable::SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const KK_FLOAT a = 1.0 - b;
    fpair = a * d_table_const.f(tidx,itable) + b * d_table_const.f(tidx,itable+1) +
      ((a*a*a-a)*d_table_const.f2(tidx,itable) + (b*b*b-b)*d_table_const.f2(tidx,itable+1)) *
      d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const KK_FLOAT fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  }
  return fpair;
}

template<ExecutionSpace Space, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
static KK_FLOAT
compute_evdwl(
    KK_FLOAT rsq,
    int itype, int jtype,
    typename PairTableRXKokkos<Space>::TableDeviceConst const& d_table_const
    ) {
  KK_FLOAT evdwl;
  Pair::union_int_float_t rsq_lookup;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (TABSTYLE == PairTable::LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    evdwl = d_table_const.e(tidx,itable);
  } else if (TABSTYLE == PairTable::LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  } else if (TABSTYLE == PairTable::SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const KK_FLOAT b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const KK_FLOAT a = 1.0 - b;
    evdwl = a * d_table_const.e(tidx,itable) + b * d_table_const.e(tidx,itable+1) +
        ((a*a*a-a)*d_table_const.e2(tidx,itable) + (b*b*b-b)*d_table_const.e2(tidx,itable+1)) *
        d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const KK_FLOAT fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  }
  return evdwl;
}

template<ExecutionSpace Space, int NEIGHFLAG, int TABSTYLE, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void
ev_tally(
    int eflag,
    int eflag_atom,
    int vflag,
    int vflag_global,
    int vflag_atom,
    int nlocal,
    int i, int j,
    EV_FLOAT& ev,
    KK_FLOAT epair, KK_FLOAT fpair,
    KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d_6::data_type,
                 typename ArrayTypes<Space>::t_float_1d_6::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_vatom,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_eatom)
{
  if (eflag) {
    if (eflag_atom) {
      auto epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) v_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) v_eatom[j] += epairhalf;
      } else {
        v_eatom[i] += epairhalf;
      }
    }
  }

  if (vflag) {
    auto v0 = delx*delx*fpair;
    auto v1 = dely*dely*fpair;
    auto v2 = delz*delz*fpair;
    auto v3 = delx*dely*fpair;
    auto v4 = delx*delz*fpair;
    auto v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR) {
          ev.v[0] += v0;
          ev.v[1] += v1;
          ev.v[2] += v2;
          ev.v[3] += v3;
          ev.v[4] += v4;
          ev.v[5] += v5;
        } else {
          if (i < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
          if (j < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
        }
      } else {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          v_vatom(i,0) += 0.5*v0;
          v_vatom(i,1) += 0.5*v1;
          v_vatom(i,2) += 0.5*v2;
          v_vatom(i,3) += 0.5*v3;
          v_vatom(i,4) += 0.5*v4;
          v_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
          v_vatom(j,0) += 0.5*v0;
          v_vatom(j,1) += 0.5*v1;
          v_vatom(j,2) += 0.5*v2;
          v_vatom(j,3) += 0.5*v3;
          v_vatom(j,4) += 0.5*v4;
          v_vatom(j,5) += 0.5*v5;
        }
      } else {
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

template <ExecutionSpace Space, int NEIGHFLAG, bool STACKPARAMS, int TABSTYLE,
          int EVFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
static EV_FLOAT
compute_item(
    int ii,
    int nlocal,
    typename ArrayTypes<Space>::t_int_1d_const const& d_ilist,
    typename ArrayTypes<Space>::t_neighbors_2d_const const& d_neighbors,
    typename ArrayTypes<Space>::t_int_1d_const const& d_numneigh,
    typename ArrayTypes<Space>::t_float_1d_3_lr_randomread const& x,
    typename ArrayTypes<Space>::t_int_1d_randomread const& type,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2,
    Few<int, 4> const& special_lj,
    Few<Few<KK_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> const& m_cutsq,
    typename ArrayTypes<Space>::t_float_2d const& d_cutsq,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d_3::data_type,
      typename ArrayTypes<Space>::t_float_1d_3::array_layout,
      typename KKDevice<typename GetDeviceType<Space>::value>::value,
      Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& f,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& uCG,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& uCGnew,
    int isite1, int isite2,
    typename PairTableRXKokkos<Space>::TableDeviceConst const& d_table_const,
    int eflag,
    int eflag_atom,
    int vflag,
    int vflag_global,
    int vflag_atom,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d_6::data_type,
                 typename ArrayTypes<Space>::t_float_1d_6::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_vatom,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_eatom) {
  EV_FLOAT ev;
  auto i = d_ilist(ii);
  auto xtmp = x(i,0);
  auto ytmp = x(i,1);
  auto ztmp = x(i,2);
  auto itype = type(i);

  auto jlist = NeighListKokkos<Space>::static_neighbors_const(i,
      d_neighbors, d_numneigh);
  auto jnum = d_numneigh(i);

  KK_FLOAT uCG_i = 0.0;
  KK_FLOAT uCGnew_i = 0.0;
  KK_FLOAT fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;

  auto mixWtSite1old_i = mixWtSite1old(i);
  auto mixWtSite2old_i = mixWtSite2old(i);
  auto mixWtSite1_i = mixWtSite1(i);
  auto mixWtSite2_i = mixWtSite2(i);

  for (int jj = 0; jj < jnum; jj++) {
    auto j = jlist(jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    auto delx = xtmp - x(j,0);
    auto dely = ytmp - x(j,1);
    auto delz = ztmp - x(j,2);
    auto rsq = delx*delx + dely*dely + delz*delz;
    auto jtype = type(j);

    if(rsq < (STACKPARAMS ? m_cutsq[itype][jtype] : d_cutsq(itype,jtype))) {
      auto mixWtSite1old_j = mixWtSite1old(j);
      auto mixWtSite2old_j = mixWtSite2old(j);
      auto mixWtSite1_j = mixWtSite1(j);
      auto mixWtSite2_j = mixWtSite2(j);

      auto fpair = factor_lj * compute_fpair<Space,TABSTYLE>(
          rsq,itype,jtype,d_table_const);

      if (isite1 == isite2) fpair *= sqrt(mixWtSite1old_i * mixWtSite2old_j);
      else fpair *= (sqrt(mixWtSite1old_i * mixWtSite2old_j) +
                     sqrt(mixWtSite2old_i * mixWtSite1old_j));

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;

      auto do_half = (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) &&
                     (NEWTON_PAIR || j < nlocal);
      if (do_half) {
        f(j,0) -= delx*fpair;
        f(j,1) -= dely*fpair;
        f(j,2) -= delz*fpair;
      }

      auto evdwl = compute_evdwl<Space,TABSTYLE>(
          rsq,itype,jtype,d_table_const);

      KK_FLOAT evdwlOld;
      if (isite1 == isite2) {
        evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwl;
        evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwl;
      } else {
        evdwlOld = (sqrt(mixWtSite1old_i*mixWtSite2old_j) +
                    sqrt(mixWtSite2old_i*mixWtSite1old_j))*evdwl;
        evdwl = (sqrt(mixWtSite1_i*mixWtSite2_j) +
                 sqrt(mixWtSite2_i*mixWtSite1_j))*evdwl;
      }
      evdwlOld *= factor_lj;
      evdwl *= factor_lj;

      uCG_i += 0.5*evdwlOld;
      if (do_half) uCG(j) += 0.5*evdwlOld;

      uCGnew_i += 0.5*evdwl;
      if (do_half) uCGnew(j) += 0.5*evdwl;
      evdwl = evdwlOld;

      ev.evdwl += (do_half ? 1.0 : 0.5)*evdwl;

      if (EVFLAG) {
        ev_tally<Space,NEIGHFLAG,TABSTYLE,NEWTON_PAIR>(
            eflag,eflag_atom,
            vflag,vflag_global,vflag_atom,
            nlocal,i,j,ev,evdwl,fpair,delx,dely,delz,
            v_vatom, v_eatom);
      }
    }
  }

  uCG(i) += uCG_i;
  uCGnew(i) += uCGnew_i;

  f(i,0) += fx_i;
  f(i,1) += fy_i;
  f(i,2) += fz_i;

  return ev;
}

template<ExecutionSpace Space, int NEIGHFLAG, bool STACKPARAMS, int TABSTYLE, int NEWTON_PAIR>
static void compute_all_items(
    EV_FLOAT& ev,
    int nlocal,
    int inum,
    typename ArrayTypes<Space>::t_int_1d_const d_ilist,
    typename ArrayTypes<Space>::t_neighbors_2d_const d_neighbors,
    typename ArrayTypes<Space>::t_int_1d_const d_numneigh,
    typename ArrayTypes<Space>::t_float_1d_3_lr_randomread x,
    typename ArrayTypes<Space>::t_int_1d_randomread type,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2,
    Few<int, 4> special_lj,
    Few<Few<KK_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> m_cutsq,
    typename ArrayTypes<Space>::t_float_2d d_cutsq,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d_3::data_type,
      typename ArrayTypes<Space>::t_float_1d_3::array_layout,
      typename KKDevice<typename GetDeviceType<Space>::value>::value,
      Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > f,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCG,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCGnew,
    int isite1, int isite2,
    typename PairTableRXKokkos<Space>::TableDeviceConst d_table_const,
    int eflag,
    int eflag_atom,
    int vflag,
    int vflag_global,
    int vflag_atom,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d_6::data_type,
                 typename ArrayTypes<Space>::t_float_1d_6::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom,
    Kokkos::View<typename ArrayTypes<Space>::t_float_1d::data_type,
                 typename ArrayTypes<Space>::t_float_1d::array_layout,
                 typename KKDevice<typename GetDeviceType<Space>::value>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom) {
  if (eflag || vflag) {
    Kokkos::parallel_reduce(inum,
    LAMMPS_LAMBDA(int i, EV_FLOAT& energy_virial) {
        energy_virial +=
          compute_item<Space,NEIGHFLAG,STACKPARAMS,TABSTYLE,1,NEWTON_PAIR>(
            i, nlocal, d_ilist, d_neighbors, d_numneigh, x, type,
            mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
            special_lj, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
            d_table_const, eflag, eflag_atom,
            vflag, vflag_global, vflag_atom, v_vatom, v_eatom);
    }, ev);
  } else {
    Kokkos::parallel_for(inum,
    LAMMPS_LAMBDA(int i) {
        compute_item<Space,NEIGHFLAG,STACKPARAMS,TABSTYLE,0,NEWTON_PAIR>(
            i, nlocal, d_ilist, d_neighbors, d_numneigh, x, type,
            mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
            special_lj, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
            d_table_const, eflag, eflag_atom,
            vflag, vflag_global, vflag_atom, v_vatom, v_eatom);
    });
  }
}

template<ExecutionSpace Space>
static void getAllMixingWeights(
    int ntotal,
    typename ArrayTypes<Space>::t_float_2d_randomread dvector,
    int nspecies,
    int isite1, int isite2,
    bool fractionalWeighting,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2old,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite1,
    typename ArrayTypes<Space>::t_float_1d const& mixWtSite2) {
  Kokkos::parallel_for(ntotal,
  LAMMPS_LAMBDA(int i) {
      getMixingWeights<Space>(dvector,nspecies,isite1,isite2,fractionalWeighting,
        i, mixWtSite1old(i), mixWtSite2old(i), mixWtSite1(i), mixWtSite2(i));
  });
}

template<ExecutionSpace Space>
template<int TABSTYLE>
void PairTableRXKokkos<Space>::compute_style(int eflag_in, int vflag_in)
{
  auto eflag = eflag_in;
  auto vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = DualViewHelper<Space>::view(k_eatom);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = DualViewHelper<Space>::view(k_vatom);
  }

  atomKK->sync(Space,datamask_read);
  if (eflag || vflag) atomKK->modified(Space,datamask_modify);
  else atomKK->modified(Space,F_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  auto type = DualViewHelper<Space>::view(atomKK->k_type);
  auto uCG = DualViewHelper<Space>::view(atomKK->k_uCG);
  auto uCGnew = DualViewHelper<Space>::view(atomKK->k_uCGnew);
  auto nlocal = atom->nlocal;
  Few<int, 4> special_lj_local;
  special_lj_local[0] = force->special_lj[0];
  special_lj_local[1] = force->special_lj[1];
  special_lj_local[2] = force->special_lj[2];
  special_lj_local[3] = force->special_lj[3];
  auto newton_pair = force->newton_pair;
  d_cutsq = d_table->cutsq;
  // loop over neighbors of my atoms

  const int ntotal = atom->nlocal + atom->nghost;
  if (ntotal > mixWtSite1.extent(0)) {
    mixWtSite1old = typename AT::t_float_1d("PairTableRXKokkos::mixWtSite1old", ntotal);
    mixWtSite2old = typename AT::t_float_1d("PairTableRXKokkos::mixWtSite2old", ntotal);
    mixWtSite1 = typename AT::t_float_1d("PairTableRXKokkos::mixWtSite1", ntotal);
    mixWtSite2 = typename AT::t_float_1d("PairTableRXKokkos::mixWtSite2", ntotal);
  }

  getAllMixingWeights<Space>(ntotal, DualViewHelper<Space>::view(atomKK->k_dvector),
      nspecies, isite1, isite2, fractionalWeighting,
      mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2);

  NeighListKokkos<Space>* l =
    dynamic_cast<NeighListKokkos<Space>*>(list);

  EV_FLOAT ev;
  if(atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        compute_all_items<Space,HALFTHREAD,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,HALFTHREAD,false,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == HALF) {
      if (newton_pair) {
        compute_all_items<Space,HALF,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,HALF,false,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        compute_all_items<Space,FULL,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,FULL,false,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    }
  } else {
    if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        compute_all_items<Space,HALFTHREAD,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,HALFTHREAD,true,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == HALF) {
      if (newton_pair) {
        compute_all_items<Space,HALF,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,HALF,true,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        compute_all_items<Space,FULL,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<Space,FULL,true,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    }
  }

  if (eflag) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

  if (eflag_atom) {
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }
}

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(k_table->k_nshiftbits,ntables,"Table::nshiftbits");
  memoryKK->create_kokkos(k_table->k_nmask,ntables,"Table::nmask");
  memoryKK->create_kokkos(k_table->k_innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(k_table->k_invdelta,ntables,"Table::invdelta");
  memoryKK->create_kokkos(k_table->k_deltasq6,ntables,"Table::deltasq6");

  if(tabstyle == LOOKUP) {
    memoryKK->create_kokkos(k_table->k_e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tlm1,"Table::f");
  }

  if(tabstyle == LINEAR) {
    memoryKK->create_kokkos(k_table->k_rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(k_table->k_de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(k_table->k_df,ntables,tlm1,"Table::df");
  }

  if(tabstyle == SPLINE) {
    memoryKK->create_kokkos(k_table->k_rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(k_table->k_e2,ntables,tablength,"Table::e2");
    memoryKK->create_kokkos(k_table->k_f2,ntables,tablength,"Table::f2");
  }

  if(tabstyle == BITMAP) {
    int ntable = 1 << tablength;
    memoryKK->create_kokkos(k_table->k_rsq,ntables,ntable,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,ntable,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,ntable,"Table::f");
    memoryKK->create_kokkos(k_table->k_de,ntables,ntable,"Table::de");
    memoryKK->create_kokkos(k_table->k_df,ntables,ntable,"Table::df");
    memoryKK->create_kokkos(k_table->k_drsq,ntables,ntable,"Table::drsq");
  }

  h_table->copy(k_table);
  d_table->copy(k_table);

  for(int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->nshiftbits[i] = tb->nshiftbits;
    h_table->nmask[i] = tb->nmask;
    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;
    h_table->deltasq6[i] = tb->deltasq6;

    for(int j = 0; j<h_table->rsq.extent(1); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for(int j = 0; j<h_table->drsq.extent(1); j++)
      h_table->drsq(i,j) = tb->drsq[j];
    for(int j = 0; j<h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for(int j = 0; j<h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
    for(int j = 0; j<h_table->f.extent(1); j++)
      h_table->f(i,j) = tb->f[j];
    for(int j = 0; j<h_table->df.extent(1); j++)
      h_table->df(i,j) = tb->df[j];
    for(int j = 0; j<h_table->e2.extent(1); j++)
      h_table->e2(i,j) = tb->e2[j];
    for(int j = 0; j<h_table->f2.extent(1); j++)
      h_table->f2(i,j) = tb->f2[j];
  }


  k_table->k_nshiftbits.modify_host();
  k_table->k_nshiftbits.sync_device();
  d_table_const.nshiftbits = d_table->nshiftbits;
  k_table->k_nmask.modify_host();
  k_table->k_nmask.sync_device();
  d_table_const.nmask = d_table->nmask;
  k_table->k_innersq.modify_host();
  k_table->k_innersq.sync_device();
  d_table_const.innersq = d_table->innersq;
  k_table->k_invdelta.modify_host();
  k_table->k_invdelta.sync_device();
  d_table_const.invdelta = d_table->invdelta;
  k_table->k_deltasq6.modify_host();
  k_table->k_deltasq6.sync_device();
  d_table_const.deltasq6 = d_table->deltasq6;

  if(tabstyle == LOOKUP) {
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
  }

  if(tabstyle == LINEAR) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_de.modify_host();
    k_table->k_de.sync_device();
    d_table_const.de = d_table->de;
    k_table->k_df.modify_host();
    k_table->k_df.sync_device();
    d_table_const.df = d_table->df;
  }

  if(tabstyle == SPLINE) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_e2.modify_host();
    k_table->k_e2.sync_device();
    d_table_const.e2 = d_table->e2;
    k_table->k_f2.modify_host();
    k_table->k_f2.sync_device();
    d_table_const.f2 = d_table->f2;
  }

  if(tabstyle == BITMAP) {
    k_table->k_rsq.modify_host();
    k_table->k_rsq.sync_device();
    d_table_const.rsq = d_table->rsq;
    k_table->k_e.modify_host();
    k_table->k_e.sync_device();
    d_table_const.e = d_table->e;
    k_table->k_f.modify_host();
    k_table->k_f.sync_device();
    d_table_const.f = d_table->f;
    k_table->k_de.modify_host();
    k_table->k_de.sync_device();
    d_table_const.de = d_table->de;
    k_table->k_df.modify_host();
    k_table->k_df.sync_device();
    d_table_const.df = d_table->df;
    k_table->k_drsq.modify_host();
    k_table->k_drsq.sync_device();
    d_table_const.drsq = d_table->drsq;
  }

  k_table->k_cutsq.modify_host();
  k_table->k_cutsq.sync_device();
  d_table_const.cutsq = d_table->cutsq;
  k_table->k_tabindex.modify_host();
  k_table->k_tabindex.sync_device();
  d_table_const.tabindex = d_table->tabindex;

  update_table = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memoryKK->create_kokkos(k_table->k_cutsq,cutsq,nt,nt,"pair:cutsq");
  memoryKK->create_kokkos(k_table->k_tabindex,tabindex,nt,nt,"pair:tabindex");
  h_table->copy(k_table);
  d_table->copy(k_table);
  d_table_const.cutsq = d_table->cutsq;
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(KK_FLOAT));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords
  // assert the tabulation is compatible with a specific long-range solver

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ewald") == 0) ewaldflag = 1;
    else if (strcmp(arg[iarg],"pppm") == 0) pppmflag = 1;
    else if (strcmp(arg[iarg],"msm") == 0) msmflag = 1;
    else if (strcmp(arg[iarg],"dispersion") == 0) dispersionflag = 1;
    else if (strcmp(arg[iarg],"tip4p") == 0) tip4pflag = 1;
    else if (strcmp(arg[iarg],"fractional") == 0) fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0) fractionalWeighting = false;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
  ntables = 0;
  tables = NULL;

  if (allocated) {
    memory->destroy(setflag);

    d_table_const.tabindex = d_table->tabindex = typename AT::t_int_2d();
    h_table->tabindex = HAT::t_int_2d();

    d_table_const.cutsq = d_table->cutsq = typename AT::t_float_2d();
    h_table->cutsq = HAT::t_float_2d();
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^rx")) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairTableRX requires a fix rx command.");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  nspecies = atom->nspecies_dpd;
  if(nspecies==0) error->all(FLERR,"There are no rx species specified.");
  int n;
  n = strlen(arg[4]) + 1;
  site1 = new char[n];
  strcpy(site1,arg[4]);

  int ispecies;
  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site1,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site1,"1fluid") != 0)
    error->all(FLERR,"Site1 name not recognized in pair coefficients");

  n = strlen(arg[5]) + 1;
  site2 = new char[n];
  strcpy(site2,arg[5]);

  for (ispecies = 0; ispecies < nspecies; ispecies++){
    if (strcmp(site2,&atom->dname[ispecies][0]) == 0) break;
  }
  if (ispecies == nspecies && strcmp(site2,"1fluid") != 0)
    error->all(FLERR,"Site2 name not recognized in pair coefficients");

  // set table cutoff

  if (narg == 7) tb->cut = force->numeric(FLERR,arg[6]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table
  // for BITMAP tables, file values can be in non-ascending order

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  KK_FLOAT rlo,rhi;
  if (tb->rflag == 0) {
    rlo = tb->rfile[0];
    rhi = tb->rfile[tb->ninput-1];
  } else {
    rlo = tb->rlo;
    rhi = tb->rhi;
  }
  if (tb->cut <= rlo || tb->cut > rhi)
    error->all(FLERR,"Invalid pair table cutoff");
  if (rlo <= 0.0) error->all(FLERR,"Invalid pair table cutoff");

  // match = 1 if don't need to spline read-in tables
  // this is only the case if r values needed by final tables
  //   exactly match r values read from file
  // for tabstyle SPLINE, always need to build spline tables

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ && tb->rhi == tb->cut) tb->match = 1;
  if (tabstyle == BITMAP && tb->ninput == 1 << tablength &&
      tb->rflag == BMP && tb->rhi == tb->cut) tb->match = 1;
  if (tb->rflag == BMP && tb->match == 0)
    error->all(FLERR,"Bitmapped table in file does not match requested table");

  // spline read-in values and compute r,e,f vectors within table

  if (tb->match == 0) spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      tabindex[i][j] = ntables;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Illegal pair_coeff command");
  ntables++;

  {
     if ( strcmp(site1,"1fluid") == 0 )
       isite1 = OneFluidValue;
     else {
       isite1 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site1, atom->dname[k]) == 0){
           isite1 = k;
           break;
         }
       }

       if (isite1 == nspecies) error->all(FLERR,"isite1 == nspecies");
     }

     if ( strcmp(site2,"1fluid") == 0 )
       isite2 = OneFluidValue;
     else {
       isite2 = nspecies;

       for (int k = 0; k < nspecies; k++){
         if (strcmp(site2, atom->dname[k]) == 0){
           isite2 = ispecies;
           break;
         }
       }

       if (isite2 == nspecies)
         error->all(FLERR,"isite2 == nspecies");
     }
  }

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairTableRXKokkos<Space>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if(i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
double PairTableRXKokkos<Space>::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  int itable;
  KK_FLOAT fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];

  typedef typename GetFloatType<Host>::type HOST_FLOAT;
  HOST_FLOAT mixWtSite1_i, mixWtSite1_j;
  HOST_FLOAT mixWtSite2_i, mixWtSite2_j;
  HOST_FLOAT mixWtSite1old_i, mixWtSite1old_j;
  HOST_FLOAT mixWtSite2old_i, mixWtSite2old_j;

  fraction = 0.0;
  a = 0.0;
  b = 0.0;

  atomKK->k_dvector.sync_host();
  HAT::t_float_2d_randomread h_dvector =
    atomKK->k_dvector.h_view;
  getMixingWeights<Host>(h_dvector,
      nspecies, isite1, isite2, fractionalWeighting,
      i,mixWtSite1old_i,mixWtSite2old_i,
      mixWtSite1_i,mixWtSite2_i);
  getMixingWeights<Host>(h_dvector,
      nspecies, isite1, isite2, fractionalWeighting,
      j,mixWtSite1old_j,mixWtSite2old_j,
      mixWtSite1_j,mixWtSite2_j);

  if (rsq < tb->innersq) error->one(FLERR,"Pair distance < table inner cutoff");

  if (tabstyle == LOOKUP) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fforce = factor_lj * tb->f[itable];
  } else if (tabstyle == LINEAR) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    b = (rsq - tb->rsq[itable]) * tb->invdelta;
    a = 1.0 - b;
    value = a * tb->f[itable] + b * tb->f[itable+1] +
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
      tb->deltasq6;
    fforce = factor_lj * value;
  } else {
    Pair::union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    itable = rsq_lookup.i & tb->nmask;
    itable >>= tb->nshiftbits;
    fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  }

  if (isite1 == isite2) fforce = sqrt(mixWtSite1_i*mixWtSite2_j)*fforce;
  else fforce = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*fforce;

  if (tabstyle == LOOKUP)
    phi = tb->e[itable];
  else if (tabstyle == LINEAR || tabstyle == BITMAP)
    phi = tb->e[itable] + fraction*tb->de[itable];
  else
    phi = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * tb->deltasq6;

  if (isite1 == isite2) phi = sqrt(mixWtSite1_i*mixWtSite2_j)*phi;
  else phi = (sqrt(mixWtSite1_i*mixWtSite2_j) + sqrt(mixWtSite2_i*mixWtSite1_j))*phi;

  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::compute_table(Table *tb)
{
  update_table = 1;
  PairTable::compute_table(tb);
}

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::init_style()
{
  neighbor->request(this,instance_me);
  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with lj/cut/kk");
  }
}

template<ExecutionSpace Space>
void PairTableRXKokkos<Space>::cleanup_copy() {
  // WHY needed: this prevents parent copy from deallocating any arrays
  allocated = 0;
  cutsq = NULL;
  eatom = NULL;
  vatom = NULL;
  h_table=NULL; d_table=NULL;
}

namespace LAMMPS_NS {
template class PairTableRXKokkos<Device>;
template class PairTableRXKokkos<Host>;

}

