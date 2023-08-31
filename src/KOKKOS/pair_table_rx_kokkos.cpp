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
   Contributing author: Dan Ibanez (SNL)
------------------------------------------------------------------------- */

#include "pair_table_rx_kokkos.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "kokkos.h"
#include "kokkos.h"
#include "kokkos_few.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cassert>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define OneFluidValue (-1)
#define isOneFluid(_site_) ( (_site_) == OneFluidValue )

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void getMixingWeights(
    typename ArrayTypes<DeviceType>::t_float_2d_randomread dvector,
    int nspecies,
    int isite1, int isite2,
    bool fractionalWeighting,
    int id,
    double &mixWtSite1old, double &mixWtSite2old,
    double &mixWtSite1, double &mixWtSite2) {
  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalOld;

  nTotal = 0.0;
  nTotalOld = 0.0;
  assert(id >= 0);
  assert(id < (int)dvector.extent(1));
  for (int ispecies = 0; ispecies < nspecies; ++ispecies) {
    assert(ispecies < (int)dvector.extent(0));
    nTotal += dvector(ispecies,id);
    assert(ispecies+nspecies < (int)dvector.extent(0));
    nTotalOld += dvector(ispecies+nspecies,id);
  }

  assert(isite1 >= 0);
  assert(isite1 < nspecies);
  assert(isite2 >= 0);
  assert(isite2 < nspecies);
  if (isOneFluid(isite1) == false) {
    nMoleculesOld1 = dvector(isite1+nspecies,id);
    nMolecules1 = dvector(isite1,id);
    fractionOld1 = nMoleculesOld1/nTotalOld;
    fraction1 = nMolecules1/nTotal;
  }
  if (isOneFluid(isite2) == false) {
    nMoleculesOld2 = dvector(isite2+nspecies,id);
    nMolecules2 = dvector(isite2,id);
    fractionOld2 = nMoleculesOld2/nTotalOld;
    fraction2 = nMolecules2/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)) {
    nMoleculesOFAold  = 0.0;
    nMoleculesOFA  = 0.0;
    fractionOFAold  = 0.0;
    fractionOFA  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++) {
      if (isite1 == ispecies || isite2 == ispecies) continue;
      nMoleculesOFAold += dvector(ispecies+nspecies,id);
      nMoleculesOFA += dvector(ispecies,id);
      fractionOFAold += dvector(ispecies+nspecies,id)/nTotalOld;
      fractionOFA += dvector(ispecies,id)/nTotal;
    }
    if (isOneFluid(isite1)) {
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if (isOneFluid(isite2)) {
      nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld2 = fractionOFAold;
      fraction2 = fractionOFA;
    }
  }

  if (fractionalWeighting) {
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

template<class DeviceType>
PairTableRXKokkos<DeviceType>::PairTableRXKokkos(LAMMPS *lmp) : PairTable(lmp)
{
  update_table = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK |
                  DVECTOR_MASK | UCG_MASK | UCGNEW_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK;
  h_table = new TableHost();
  d_table = new TableDevice();
  fractionalWeighting = true;

  site1 = nullptr;
  site2 = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTableRXKokkos<DeviceType>::~PairTableRXKokkos()
{
  if (copymode) return;

  delete [] site1;
  delete [] site2;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memory->destroy(setflag);
    memoryKK->destroy_kokkos(d_table->cutsq,cutsq);
    memoryKK->destroy_kokkos(d_table->tabindex,tabindex);
  }

  delete h_table;
  h_table = nullptr;
  delete d_table;
  d_table = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (update_table)
    create_kokkos_tables();
  if (tabstyle == LOOKUP)
    compute_style<LOOKUP>(eflag_in,vflag_in);
  if (tabstyle == LINEAR)
    compute_style<LINEAR>(eflag_in,vflag_in);
  if (tabstyle == SPLINE)
    compute_style<SPLINE>(eflag_in,vflag_in);
  if (tabstyle == BITMAP)
    compute_style<BITMAP>(eflag_in,vflag_in);
}

KOKKOS_INLINE_FUNCTION static int sbmask(const int& j)
{
  return j >> SBBITS & 3;
}

template <class DeviceType, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
static F_FLOAT
compute_fpair(F_FLOAT rsq,
              int itype, int jtype,
              typename PairTableRXKokkos<DeviceType>::TableDeviceConst const& d_table_const
              ) {
  Pair::union_int_float_t rsq_lookup;
  double fpair;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (TABSTYLE == PairTable::LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    fpair = d_table_const.f(tidx,itable);
  } else if (TABSTYLE == PairTable::LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  } else if (TABSTYLE == PairTable::SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const double a = 1.0 - b;
    fpair = a * d_table_const.f(tidx,itable) + b * d_table_const.f(tidx,itable+1) +
      ((a*a*a-a)*d_table_const.f2(tidx,itable) + (b*b*b-b)*d_table_const.f2(tidx,itable+1)) *
      d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const double fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    fpair = d_table_const.f(tidx,itable) + fraction*d_table_const.df(tidx,itable);
  }
  return fpair;
}

template<class DeviceType, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
static F_FLOAT
compute_evdwl(
    F_FLOAT rsq,
    int itype, int jtype,
    typename PairTableRXKokkos<DeviceType>::TableDeviceConst const& d_table_const
    ) {
  double evdwl;
  Pair::union_int_float_t rsq_lookup;
  const int tidx = d_table_const.tabindex(itype,jtype);
  if (TABSTYLE == PairTable::LOOKUP) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    evdwl = d_table_const.e(tidx,itable);
  } else if (TABSTYLE == PairTable::LINEAR) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double fraction = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  } else if (TABSTYLE == PairTable::SPLINE) {
    const int itable = static_cast<int> ((rsq - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
    const double b = (rsq - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx);
    const double a = 1.0 - b;
    evdwl = a * d_table_const.e(tidx,itable) + b * d_table_const.e(tidx,itable+1) +
        ((a*a*a-a)*d_table_const.e2(tidx,itable) + (b*b*b-b)*d_table_const.e2(tidx,itable+1)) *
        d_table_const.deltasq6(tidx);
  } else {
    rsq_lookup.f = rsq;
    int itable = rsq_lookup.i & d_table_const.nmask(tidx);
    itable >>= d_table_const.nshiftbits(tidx);
    const double fraction = (rsq_lookup.f - d_table_const.rsq(tidx,itable)) * d_table_const.drsq(tidx,itable);
    evdwl = d_table_const.e(tidx,itable) + fraction*d_table_const.de(tidx,itable);
  }
  return evdwl;
}

template<class DeviceType, int NEIGHFLAG, int TABSTYLE, int NEWTON_PAIR>
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
    F_FLOAT epair, F_FLOAT fpair,
    F_FLOAT delx, F_FLOAT dely, F_FLOAT delz,
    Kokkos::View<F_FLOAT*[6],
                 typename ArrayTypes<DeviceType>::t_virial_array::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_vatom,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
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

template <class DeviceType, int NEIGHFLAG, bool STACKPARAMS, int TABSTYLE,
          int EVFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
static EV_FLOAT
compute_item(
    int ii,
    int nlocal,
    typename ArrayTypes<DeviceType>::t_int_1d_const const& d_ilist,
    typename ArrayTypes<DeviceType>::t_neighbors_2d_const const& d_neighbors,
    typename ArrayTypes<DeviceType>::t_int_1d_const const& d_numneigh,
    typename ArrayTypes<DeviceType>::t_x_array_randomread const& x,
    typename ArrayTypes<DeviceType>::t_int_1d_randomread const& type,
    Kokkos::View<double*, DeviceType> const& mixWtSite1old,
    Kokkos::View<double*, DeviceType> const& mixWtSite2old,
    Kokkos::View<double*, DeviceType> const& mixWtSite1,
    Kokkos::View<double*, DeviceType> const& mixWtSite2,
    Few<int, 4> const& special_lj,
    Few<Few<F_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> const& m_cutsq,
    typename ArrayTypes<DeviceType>::t_ffloat_2d const& d_cutsq,
    Kokkos::View<F_FLOAT*[3],
      typename ArrayTypes<DeviceType>::t_f_array::array_layout,
      typename KKDevice<DeviceType>::value,
      Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& f,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& uCG,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& uCGnew,
    int isite1, int isite2,
    typename PairTableRXKokkos<DeviceType>::TableDeviceConst const& d_table_const,
    int eflag,
    int eflag_atom,
    int vflag,
    int vflag_global,
    int vflag_atom,
    Kokkos::View<F_FLOAT*[6],
                 typename ArrayTypes<DeviceType>::t_virial_array::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_vatom,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > const& v_eatom) {
  EV_FLOAT ev;
  auto i = d_ilist(ii);
  auto xtmp = x(i,0);
  auto ytmp = x(i,1);
  auto ztmp = x(i,2);
  auto itype = type(i);

  auto jlist = NeighListKokkos<DeviceType>::static_neighbors_const(i,
      d_neighbors, d_numneigh);
  auto jnum = d_numneigh(i);

  double uCG_i = 0.0;
  double uCGnew_i = 0.0;
  double fx_i = 0.0, fy_i = 0.0, fz_i = 0.0;

  auto mixWtSite1old_i = mixWtSite1old(i);
  auto mixWtSite2old_i = mixWtSite2old(i);
  auto mixWtSite1_i = mixWtSite1(i);
  auto mixWtSite2_i = mixWtSite2(i);

  for (int jj = 0; jj < jnum; jj++) {
    auto j = jlist(jj);
    const F_FLOAT factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    auto delx = xtmp - x(j,0);
    auto dely = ytmp - x(j,1);
    auto delz = ztmp - x(j,2);
    auto rsq = delx*delx + dely*dely + delz*delz;
    auto jtype = type(j);

    if (rsq < (STACKPARAMS ? m_cutsq[itype][jtype] : d_cutsq(itype,jtype))) {
      auto mixWtSite1old_j = mixWtSite1old(j);
      auto mixWtSite2old_j = mixWtSite2old(j);
      auto mixWtSite1_j = mixWtSite1(j);
      auto mixWtSite2_j = mixWtSite2(j);

      auto fpair = factor_lj * compute_fpair<DeviceType,TABSTYLE>(
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

      auto evdwl = compute_evdwl<DeviceType,TABSTYLE>(
          rsq,itype,jtype,d_table_const);

      double evdwlOld;
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
        ev_tally<DeviceType,NEIGHFLAG,TABSTYLE,NEWTON_PAIR>(
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

template<class DeviceType, int NEIGHFLAG, bool STACKPARAMS, int TABSTYLE, int NEWTON_PAIR>
static void compute_all_items(
    EV_FLOAT& ev,
    int nlocal,
    int inum,
    typename ArrayTypes<DeviceType>::t_int_1d_const d_ilist,
    typename ArrayTypes<DeviceType>::t_neighbors_2d_const d_neighbors,
    typename ArrayTypes<DeviceType>::t_int_1d_const d_numneigh,
    typename ArrayTypes<DeviceType>::t_x_array_randomread x,
    typename ArrayTypes<DeviceType>::t_int_1d_randomread type,
    Kokkos::View<double*, DeviceType> const& mixWtSite1old,
    Kokkos::View<double*, DeviceType> const& mixWtSite2old,
    Kokkos::View<double*, DeviceType> const& mixWtSite1,
    Kokkos::View<double*, DeviceType> const& mixWtSite2,
    Few<int, 4> special_lj,
    Few<Few<F_FLOAT, MAX_TYPES_STACKPARAMS+1>, MAX_TYPES_STACKPARAMS+1> m_cutsq,
    typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq,
    Kokkos::View<F_FLOAT*[3],
      typename ArrayTypes<DeviceType>::t_f_array::array_layout,
      typename KKDevice<DeviceType>::value,
      Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > f,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCG,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > uCGnew,
    int isite1, int isite2,
    typename PairTableRXKokkos<DeviceType>::TableDeviceConst d_table_const,
    int eflag,
    int eflag_atom,
    int vflag,
    int vflag_global,
    int vflag_atom,
    Kokkos::View<F_FLOAT*[6],
                 typename ArrayTypes<DeviceType>::t_virial_array::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom,
    Kokkos::View<E_FLOAT*,
                 typename ArrayTypes<DeviceType>::t_efloat_1d::array_layout,
                 typename KKDevice<DeviceType>::value,
                 Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom) {
  if (eflag || vflag) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,inum),
     LAMMPS_LAMBDA(int i, EV_FLOAT& energy_virial) {
        energy_virial +=
          compute_item<DeviceType,NEIGHFLAG,STACKPARAMS,TABSTYLE,1,NEWTON_PAIR>(
            i, nlocal, d_ilist, d_neighbors, d_numneigh, x, type,
            mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
            special_lj, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
            d_table_const, eflag, eflag_atom,
            vflag, vflag_global, vflag_atom, v_vatom, v_eatom);
    }, ev);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),
     LAMMPS_LAMBDA(int i) {
        compute_item<DeviceType,NEIGHFLAG,STACKPARAMS,TABSTYLE,0,NEWTON_PAIR>(
            i, nlocal, d_ilist, d_neighbors, d_numneigh, x, type,
            mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
            special_lj, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
            d_table_const, eflag, eflag_atom,
            vflag, vflag_global, vflag_atom, v_vatom, v_eatom);
    });
  }
}

template<class DeviceType>
static void getAllMixingWeights(
    int ntotal,
    typename ArrayTypes<DeviceType>::t_float_2d_randomread dvector,
    int nspecies,
    int isite1, int isite2,
    bool fractionalWeighting,
    Kokkos::View<double*, DeviceType> const& mixWtSite1old,
    Kokkos::View<double*, DeviceType> const& mixWtSite2old,
    Kokkos::View<double*, DeviceType> const& mixWtSite1,
    Kokkos::View<double*, DeviceType> const& mixWtSite2) {
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,ntotal),
   LAMMPS_LAMBDA(int i) {
      getMixingWeights<DeviceType>(dvector,nspecies,isite1,isite2,fractionalWeighting,
        i, mixWtSite1old(i), mixWtSite2old(i), mixWtSite1(i), mixWtSite2(i));
  });
}

template<class DeviceType>
template<int TABSTYLE>
void PairTableRXKokkos<DeviceType>::compute_style(int eflag_in, int vflag_in)
{
  auto eflag = eflag_in;
  auto vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  auto type = atomKK->k_type.view<DeviceType>();
  auto uCG = atomKK->k_uCG.view<DeviceType>();
  auto uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  auto nlocal = atom->nlocal;
  Few<int, 4> special_lj_local;
  special_lj_local[0] = force->special_lj[0];
  special_lj_local[1] = force->special_lj[1];
  special_lj_local[2] = force->special_lj[2];
  special_lj_local[3] = force->special_lj[3];
  auto newton_pair = force->newton_pair;
  d_cutsq = d_table->cutsq;
  // loop over neighbors of my atoms

  copymode = 1;

  const int ntotal = atom->nlocal + atom->nghost;
  if (ntotal > (int)mixWtSite1.extent(0)) {
    mixWtSite1old = Kokkos::View<double*, DeviceType>("PairTableRXKokkos::mixWtSite1old", ntotal);
    mixWtSite2old = Kokkos::View<double*, DeviceType>("PairTableRXKokkos::mixWtSite2old", ntotal);
    mixWtSite1 = Kokkos::View<double*, DeviceType>("PairTableRXKokkos::mixWtSite1", ntotal);
    mixWtSite2 = Kokkos::View<double*, DeviceType>("PairTableRXKokkos::mixWtSite2", ntotal);
  }

  getAllMixingWeights(ntotal, atomKK->k_dvector.template view<DeviceType>(),
      nspecies, isite1, isite2, fractionalWeighting,
      mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2);

  NeighListKokkos<DeviceType>* l =
    dynamic_cast<NeighListKokkos<DeviceType>*>(list);

  EV_FLOAT ev;
  if (atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        compute_all_items<DeviceType,HALFTHREAD,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,HALFTHREAD,false,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == HALF) {
      if (newton_pair) {
        compute_all_items<DeviceType,HALF,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,HALF,false,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        compute_all_items<DeviceType,FULL,false,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,FULL,false,TABSTYLE,0>(
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
        compute_all_items<DeviceType,HALFTHREAD,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,HALFTHREAD,true,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == HALF) {
      if (newton_pair) {
        compute_all_items<DeviceType,HALF,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,HALF,true,TABSTYLE,0>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        compute_all_items<DeviceType,FULL,true,TABSTYLE,1>(
          ev, nlocal, l->inum, l->d_ilist, l->d_neighbors, l->d_numneigh,
          x, type, mixWtSite1old, mixWtSite2old, mixWtSite1, mixWtSite2,
          special_lj_local, m_cutsq, d_cutsq, f, uCG, uCGnew, isite1, isite2,
          d_table_const, eflag, eflag_atom,
          vflag, vflag_global, vflag_atom, d_vatom, d_eatom);
      } else {
        compute_all_items<DeviceType,FULL,true,TABSTYLE,0>(
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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;
}

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(d_table->nshiftbits,h_table->nshiftbits,ntables,"Table::nshiftbits");
  memoryKK->create_kokkos(d_table->nmask,h_table->nmask,ntables,"Table::nmask");
  memoryKK->create_kokkos(d_table->innersq,h_table->innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(d_table->invdelta,h_table->invdelta,ntables,"Table::invdelta");
  memoryKK->create_kokkos(d_table->deltasq6,h_table->deltasq6,ntables,"Table::deltasq6");

  if (tabstyle == LOOKUP) {
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tlm1,"Table::f");
  }

  if (tabstyle == LINEAR) {
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(d_table->df,h_table->df,ntables,tlm1,"Table::df");
  }

  if (tabstyle == SPLINE) {
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(d_table->e2,h_table->e2,ntables,tablength,"Table::e2");
    memoryKK->create_kokkos(d_table->f2,h_table->f2,ntables,tablength,"Table::f2");
  }

  if (tabstyle == BITMAP) {
    int ntable = 1 << tablength;
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,ntable,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,ntable,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,ntable,"Table::f");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,ntable,"Table::de");
    memoryKK->create_kokkos(d_table->df,h_table->df,ntables,ntable,"Table::df");
    memoryKK->create_kokkos(d_table->drsq,h_table->drsq,ntables,ntable,"Table::drsq");
  }



  for (int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->nshiftbits[i] = tb->nshiftbits;
    h_table->nmask[i] = tb->nmask;
    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;
    h_table->deltasq6[i] = tb->deltasq6;

    for (int j = 0; j < (int)h_table->rsq.extent(1); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for (int j = 0; j < (int)h_table->drsq.extent(1); j++)
      h_table->drsq(i,j) = tb->drsq[j];
    for (int j = 0; j < (int)h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for (int j = 0; j < (int)h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
    for (int j = 0; j < (int)h_table->f.extent(1); j++)
      h_table->f(i,j) = tb->f[j];
    for (int j = 0; j < (int)h_table->df.extent(1); j++)
      h_table->df(i,j) = tb->df[j];
    for (int j = 0; j < (int)h_table->e2.extent(1); j++)
      h_table->e2(i,j) = tb->e2[j];
    for (int j = 0; j < (int)h_table->f2.extent(1); j++)
      h_table->f2(i,j) = tb->f2[j];
  }


  Kokkos::deep_copy(d_table->nshiftbits,h_table->nshiftbits);
  d_table_const.nshiftbits = d_table->nshiftbits;
  Kokkos::deep_copy(d_table->nmask,h_table->nmask);
  d_table_const.nmask = d_table->nmask;
  Kokkos::deep_copy(d_table->innersq,h_table->innersq);
  d_table_const.innersq = d_table->innersq;
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  d_table_const.invdelta = d_table->invdelta;
  Kokkos::deep_copy(d_table->deltasq6,h_table->deltasq6);
  d_table_const.deltasq6 = d_table->deltasq6;

  if (tabstyle == LOOKUP) {
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
  }

  if (tabstyle == LINEAR) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->de,h_table->de);
    d_table_const.de = d_table->de;
    Kokkos::deep_copy(d_table->df,h_table->df);
    d_table_const.df = d_table->df;
  }

  if (tabstyle == SPLINE) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->e2,h_table->e2);
    d_table_const.e2 = d_table->e2;
    Kokkos::deep_copy(d_table->f2,h_table->f2);
    d_table_const.f2 = d_table->f2;
  }

  if (tabstyle == BITMAP) {
    Kokkos::deep_copy(d_table->rsq,h_table->rsq);
    d_table_const.rsq = d_table->rsq;
    Kokkos::deep_copy(d_table->e,h_table->e);
    d_table_const.e = d_table->e;
    Kokkos::deep_copy(d_table->f,h_table->f);
    d_table_const.f = d_table->f;
    Kokkos::deep_copy(d_table->de,h_table->de);
    d_table_const.de = d_table->de;
    Kokkos::deep_copy(d_table->df,h_table->df);
    d_table_const.df = d_table->df;
    Kokkos::deep_copy(d_table->drsq,h_table->drsq);
    d_table_const.drsq = d_table->drsq;
  }

  Kokkos::deep_copy(d_table->cutsq,h_table->cutsq);
  d_table_const.cutsq = d_table->cutsq;
  Kokkos::deep_copy(d_table->tabindex,h_table->tabindex);
  d_table_const.tabindex = d_table->tabindex;

  update_table = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memoryKK->create_kokkos(d_table->cutsq,h_table->cutsq,cutsq,nt,nt,"pair:cutsq");
  memoryKK->create_kokkos(d_table->tabindex,h_table->tabindex,tabindex,nt,nt,"pair:tabindex");
  d_table_const.cutsq = d_table->cutsq;
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(double));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = utils::inumeric(FLERR,arg[1],false,lmp);
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
  tables = nullptr;

  if (allocated) {
    memory->destroy(setflag);

    d_table_const.tabindex = d_table->tabindex = typename ArrayTypes<DeviceType>::t_int_2d();
    h_table->tabindex = typename ArrayTypes<LMPHostType>::t_int_2d();

    d_table_const.cutsq = d_table->cutsq = typename ArrayTypes<DeviceType>::t_ffloat_2d();
    h_table->cutsq = typename ArrayTypes<LMPHostType>::t_ffloat_2d();
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::coeff(int narg, char **arg)
{
  if (narg != 6 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  bool rx_flag = false;
  for (int i = 0; i < modify->nfix; i++)
    if (utils::strmatch(modify->fix[i]->style,"^rx")) rx_flag = true;
  if (!rx_flag) error->all(FLERR,"PairTableRX requires a fix rx command.");

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  nspecies = atom->nspecies_dpd;
  if (nspecies==0) error->all(FLERR,"There are no rx species specified.");
  site1 = utils::strdup(arg[4]);

  int ispecies;
  for (ispecies = 0; ispecies < nspecies; ispecies++) {
    if (strcmp(site1,&atom->dvname[ispecies][0]) == 0) break;
  }

  if (ispecies == nspecies && strcmp(site1,"1fluid") != 0)
    error->all(FLERR,"Site1 name not recognized in pair coefficients");
  site2 = utils::strdup(arg[5]);

  for (ispecies = 0; ispecies < nspecies; ispecies++)
    if (strcmp(site2,&atom->dvname[ispecies][0]) == 0) break;

  if (ispecies == nspecies && strcmp(site2,"1fluid") != 0)
    error->all(FLERR,"Site2 name not recognized in pair coefficients");

  // set table cutoff

  if (narg == 7) tb->cut = utils::numeric(FLERR,arg[6],false,lmp);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // ensure cutoff is within table
  // for BITMAP tables, file values can be in non-ascending order

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  double rlo,rhi;
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
     if (strcmp(site1,"1fluid") == 0)
       isite1 = OneFluidValue;
     else {
       isite1 = nspecies;

       for (int k = 0; k < nspecies; k++) {
         if (strcmp(site1, atom->dvname[k]) == 0) {
           isite1 = k;
           break;
         }
       }

       if (isite1 == nspecies) error->all(FLERR,"isite1 == nspecies");
     }

     if (strcmp(site2,"1fluid") == 0)
       isite2 = OneFluidValue;
     else {
       isite2 = nspecies;

       for (int k = 0; k < nspecies; k++) {
         if (strcmp(site2, atom->dvname[k]) == 0) {
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

template<class DeviceType>
double PairTableRXKokkos<DeviceType>::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_cutsq[j][i] = m_cutsq[i][j] = tables[tabindex[i][j]].cut*tables[tabindex[i][j]].cut;
  }

  return tables[tabindex[i][j]].cut;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairTableRXKokkos<DeviceType>::single(int i, int j, int itype, int jtype, double rsq,
                                             double /*factor_coul*/, double factor_lj,
                         double &fforce)
{
  int itable;
  double fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;
  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;

  fraction = 0.0;
  a = 0.0;
  b = 0.0;

  atomKK->k_dvector.template sync<LMPHostType>();
  typename ArrayTypes<LMPHostType>::t_float_2d_randomread h_dvector =
    atomKK->k_dvector.view<LMPHostType>();
  getMixingWeights<LMPHostType>(h_dvector,
      nspecies, isite1, isite2, fractionalWeighting,
      i,mixWtSite1old_i,mixWtSite2old_i,
      mixWtSite1_i,mixWtSite2_i);
  getMixingWeights<LMPHostType>(h_dvector,
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

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::compute_table(Table *tb)
{
  update_table = 1;
  PairTable::compute_table(tb);
}

template<class DeviceType>
void PairTableRXKokkos<DeviceType>::init_style()
{
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->add_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

namespace LAMMPS_NS {
template class PairTableRXKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairTableRXKokkos<LMPHostType>;
#endif

}

