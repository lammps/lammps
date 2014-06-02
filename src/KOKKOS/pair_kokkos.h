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

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_KOKKOS_H
#define LMP_PAIR_KOKKOS_H

#include "Kokkos_Macros.hpp"
#include "pair.h"
#include "neigh_list_kokkos.h"
#include "Kokkos_Vectorization.hpp"

namespace LAMMPS_NS {

template <class PairStyle, int NEIGHFLAG, bool STACKPARAMS, class Specialisation = void>
struct PairComputeFunctor  {
  typedef typename PairStyle::device_type device_type ;
  typedef EV_FLOAT value_type;

  PairStyle c;
  NeighListKokkos<device_type> list;

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr) {};
  ~PairComputeFunctor() {c.cleanup_copy();list.clean_copy();};

  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list) const {
    EV_FLOAT ev;
    const int i = list.d_ilist[ii];
    const X_FLOAT xtmp = c.x(i,0);
    const X_FLOAT ytmp = c.x(i,1);
    const X_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = neighbors_i(jj);
      const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const X_FLOAT delx = xtmp - c.x(j,0);
      const X_FLOAT dely = ytmp - c.x(j,1);
      const X_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if ((NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < c.nlocal)) {
          Kokkos::atomic_fetch_add(&c.f(j,0),-delx*fpair);
          Kokkos::atomic_fetch_add(&c.f(j,1),-dely*fpair);
          Kokkos::atomic_fetch_add(&c.f(j,2),-delz*fpair);
        }

        if ((NEIGHFLAG==HALF) && (NEWTON_PAIR || j < c.nlocal)) {
          c.f(j,0) -= delx*fpair;
          c.f(j,1) -= dely*fpair;
          c.f(j,2) -= delz*fpair;
        }

        if (EVFLAG) {
          if (c.eflag) {
            ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?1.0:0.5)*
              factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            if (c.COUL_FLAG)
              ev.ecoul += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<c.nlocal)))?1.0:0.5)*
                factor_lj * c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          }

          if (c.vflag_either) ev_tally(ev,i,j,fpair,delx,dely,delz);
        }
      }

    }
    if (NEIGHFLAG == HALFTHREAD) {
      Kokkos::atomic_fetch_add(&c.f(i,0),fxtmp);
      Kokkos::atomic_fetch_add(&c.f(i,1),fytmp);
      Kokkos::atomic_fetch_add(&c.f(i,2),fztmp);
    } else {
      c.f(i,0) += fxtmp;
      c.f(i,1) += fytmp;
      c.f(i,2) += fztmp;
    }

    return ev;
  }

  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const
  {
    const int EFLAG = c.eflag;
    const int NEWTON_PAIR = c.newton_pair;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const E_FLOAT epairhalf = 0.5 * (ev.evdwl + ev.ecoul);
        if (NEWTON_PAIR || i < c.nlocal) c.eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < c.nlocal) c.eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const E_FLOAT v0 = delx*delx*fpair;
      const E_FLOAT v1 = dely*dely*fpair;
      const E_FLOAT v2 = delz*delz*fpair;
      const E_FLOAT v3 = delx*dely*fpair;
      const E_FLOAT v4 = delx*delz*fpair;
      const E_FLOAT v5 = dely*delz*fpair;

      if (c.vflag_global) {
        if (NEIGHFLAG) {
          if (NEWTON_PAIR) {
            ev.v[0] += v0;
            ev.v[1] += v1;
            ev.v[2] += v2;
            ev.v[3] += v3;
            ev.v[4] += v4;
            ev.v[5] += v5;
          } else {
            if (i < c.nlocal) {
              ev.v[0] += 0.5*v0;
              ev.v[1] += 0.5*v1;
              ev.v[2] += 0.5*v2;
              ev.v[3] += 0.5*v3;
              ev.v[4] += 0.5*v4;
              ev.v[5] += 0.5*v5;
            }
            if (j < c.nlocal) {
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

      if (c.vflag_atom) {
        if (NEWTON_PAIR || i < c.nlocal) {
          c.d_vatom(i,0) += 0.5*v0;
          c.d_vatom(i,1) += 0.5*v1;
          c.d_vatom(i,2) += 0.5*v2;
          c.d_vatom(i,3) += 0.5*v3;
          c.d_vatom(i,4) += 0.5*v4;
          c.d_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || (NEIGHFLAG && j < c.nlocal)) {
        c.d_vatom(j,0) += 0.5*v0;
        c.d_vatom(j,1) += 0.5*v1;
        c.d_vatom(j,2) += 0.5*v2;
        c.d_vatom(j,3) += 0.5*v3;
        c.d_vatom(j,4) += 0.5*v4;
        c.d_vatom(j,5) += 0.5*v5;
        }
      }
    }
  }


  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (c.newton_pair) compute_item<0,1>(i,list);
    else compute_item<0,0>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    if (c.newton_pair)
      energy_virial += compute_item<1,1>(i,list);
    else
      energy_virial += compute_item<1,0>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  static void init(volatile value_type &update) {
    update.evdwl = 0;
    update.ecoul = 0;
    update.v[0] = 0;
    update.v[1] = 0;
    update.v[2] = 0;
    update.v[3] = 0;
    update.v[4] = 0;
    update.v[5] = 0;
  }
  KOKKOS_INLINE_FUNCTION 
  static void join(volatile value_type &update,
                   const volatile value_type &source) {
    update.evdwl += source.evdwl;
    update.ecoul += source.ecoul;
    update.v[0] += source.v[0];
    update.v[1] += source.v[1];
    update.v[2] += source.v[2];
    update.v[3] += source.v[3];
    update.v[4] += source.v[4];
    update.v[5] += source.v[5];
  }


};

template <class PairStyle, bool STACKPARAMS, class Specialisation>
struct PairComputeFunctor<PairStyle,FULLCLUSTER,STACKPARAMS,Specialisation>  {
  typedef typename PairStyle::device_type device_type ;
  typedef Kokkos::Vectorization<device_type,NeighClusterSize> vectorization;
  typedef EV_FLOAT value_type;

  PairStyle c;
  NeighListKokkos<device_type> list;

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr) {};
  ~PairComputeFunctor() {c.cleanup_copy();list.clean_copy();};

  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const device_type& dev,
                        const NeighListKokkos<device_type> &list) const {
    EV_FLOAT ev;
    const int i = vectorization::global_thread_rank(dev);

    const X_FLOAT xtmp = c.c_x(i,0);
    const X_FLOAT ytmp = c.c_x(i,1);
    const X_FLOAT ztmp = c.c_x(i,2);
    const int itype = c.type(i);

    const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = list.d_numneigh[i];

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      const int jjj = neighbors_i(jj);

      for (int k = vectorization::begin(); k<NeighClusterSize; k+=vectorization::increment) {
        const F_FLOAT factor_lj = c.special_lj[sbmask(jjj+k)];
        const int j = (jjj + k)&NEIGHMASK;
        if((j==i)||(j>=c.nall)) continue;
        const X_FLOAT delx = xtmp - c.c_x(j,0);
        const X_FLOAT dely = ytmp - c.c_x(j,1);
        const X_FLOAT delz = ztmp - c.c_x(j,2);
        const int jtype = c.type(j);
        const F_FLOAT rsq = (delx*delx + dely*dely + delz*delz);

        if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

          const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          fxtmp += delx*fpair;
          fytmp += dely*fpair;
          fztmp += delz*fpair;

          if (EVFLAG) {
            if (c.eflag) {
              ev.evdwl += 0.5*
                factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
              if (c.COUL_FLAG)
                ev.ecoul += 0.5*
                  factor_lj * c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            }

            if (c.vflag_either) ev_tally(ev,i,j,fpair,delx,dely,delz);
          }
        }
      }
    }

    const F_FLOAT fx = vectorization::reduce(fxtmp);
    const F_FLOAT fy = vectorization::reduce(fytmp);
    const F_FLOAT fz = vectorization::reduce(fztmp);
    if(vectorization::is_lane_0(dev)) {
      c.f(i,0) += fx;
      c.f(i,1) += fy;
      c.f(i,2) += fz;
    }

    return ev;
  }

  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const
  {
    const int EFLAG = c.eflag;
    const int NEWTON_PAIR = c.newton_pair;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const E_FLOAT epairhalf = 0.5 * (ev.evdwl + ev.ecoul);
        if (NEWTON_PAIR || i < c.nlocal) c.eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < c.nlocal) c.eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const E_FLOAT v0 = delx*delx*fpair;
      const E_FLOAT v1 = dely*dely*fpair;
      const E_FLOAT v2 = delz*delz*fpair;
      const E_FLOAT v3 = delx*dely*fpair;
      const E_FLOAT v4 = delx*delz*fpair;
      const E_FLOAT v5 = dely*delz*fpair;

      if (c.vflag_global) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
      }

      if (c.vflag_atom) {
        if (i < c.nlocal) {
          c.d_vatom(i,0) += 0.5*v0;
          c.d_vatom(i,1) += 0.5*v1;
          c.d_vatom(i,2) += 0.5*v2;
          c.d_vatom(i,3) += 0.5*v3;
          c.d_vatom(i,4) += 0.5*v4;
          c.d_vatom(i,5) += 0.5*v5;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const device_type& dev) const {
    if (c.newton_pair) compute_item<0,1>(dev,list);
    else compute_item<0,0>(dev,list);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const device_type& dev, value_type &energy_virial) const {
    if (c.newton_pair)
      energy_virial += compute_item<1,1>(dev,list);
    else
      energy_virial += compute_item<1,0>(dev,list);
  }

  KOKKOS_INLINE_FUNCTION
  static void init(volatile value_type &update) {
    update.evdwl = 0;
    update.ecoul = 0;
    update.v[0] = 0;
    update.v[1] = 0;
    update.v[2] = 0;
    update.v[3] = 0;
    update.v[4] = 0;
    update.v[5] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &source) {
    update.evdwl += source.evdwl;
    update.ecoul += source.ecoul;
    update.v[0] += source.v[0];
    update.v[1] += source.v[1];
    update.v[2] += source.v[2];
    update.v[3] += source.v[3];
    update.v[4] += source.v[4];
    update.v[5] += source.v[5];
  }


};

template <class PairStyle, bool STACKPARAMS, class Specialisation>
struct PairComputeFunctor<PairStyle,N2,STACKPARAMS,Specialisation>  {
  typedef typename PairStyle::device_type device_type ;
  typedef EV_FLOAT value_type;

  PairStyle c;
  NeighListKokkos<device_type> list;

  PairComputeFunctor(PairStyle* c_ptr,
                          NeighListKokkos<device_type>* list_ptr):
  c(*c_ptr),list(*list_ptr) {};
  ~PairComputeFunctor() {c.cleanup_copy();list.clean_copy();};

  KOKKOS_INLINE_FUNCTION int sbmask(const int& j) const {
    return j >> SBBITS & 3;
  }

  template<int EVFLAG, int NEWTON_PAIR>
  KOKKOS_FUNCTION
  EV_FLOAT compute_item(const int& ii,
                        const NeighListKokkos<device_type> &list) const {
    EV_FLOAT ev;
    const int i = ii;//list.d_ilist[ii];
    const X_FLOAT xtmp = c.x(i,0);
    const X_FLOAT ytmp = c.x(i,1);
    const X_FLOAT ztmp = c.x(i,2);
    const int itype = c.type(i);

    //const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
    const int jnum = c.nall;

    F_FLOAT fxtmp = 0.0;
    F_FLOAT fytmp = 0.0;
    F_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jj;//neighbors_i(jj);
      if(i==j) continue;
      const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
      j &= NEIGHMASK;
      const X_FLOAT delx = xtmp - c.x(j,0);
      const X_FLOAT dely = ytmp - c.x(j,1);
      const X_FLOAT delz = ztmp - c.x(j,2);
      const int jtype = c.type(j);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if(rsq < (STACKPARAMS?c.m_cutsq[itype][jtype]:c.d_cutsq(itype,jtype))) {

        const F_FLOAT fpair = factor_lj*c.template compute_fpair<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;

        if (EVFLAG) {
          if (c.eflag) {
            ev.evdwl += 0.5*
              factor_lj * c.template compute_evdwl<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
            if (c.COUL_FLAG)
              ev.ecoul += 0.5*
                factor_lj * c.template compute_ecoul<STACKPARAMS,Specialisation>(rsq,i,j,itype,jtype);
          }

          if (c.vflag_either) ev_tally(ev,i,j,fpair,delx,dely,delz);
        }
      }
    }

    c.f(i,0) += fxtmp;
    c.f(i,1) += fytmp;
    c.f(i,2) += fztmp;

    return ev;
  }

  KOKKOS_INLINE_FUNCTION
    void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const
  {
    const int EFLAG = c.eflag;
    const int VFLAG = c.vflag_either;

    if (EFLAG) {
      if (c.eflag_atom) {
        const E_FLOAT epairhalf = 0.5 * (ev.evdwl + ev.ecoul);
        if (i < c.nlocal) c.eatom[i] += epairhalf;
        if (j < c.nlocal) c.eatom[j] += epairhalf;
      }
    }

    if (VFLAG) {
      const E_FLOAT v0 = delx*delx*fpair;
      const E_FLOAT v1 = dely*dely*fpair;
      const E_FLOAT v2 = delz*delz*fpair;
      const E_FLOAT v3 = delx*dely*fpair;
      const E_FLOAT v4 = delx*delz*fpair;
      const E_FLOAT v5 = dely*delz*fpair;

      if (c.vflag_global) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
      }

      if (c.vflag_atom) {
        if (i < c.nlocal) {
          c.d_vatom(i,0) += 0.5*v0;
          c.d_vatom(i,1) += 0.5*v1;
          c.d_vatom(i,2) += 0.5*v2;
          c.d_vatom(i,3) += 0.5*v3;
          c.d_vatom(i,4) += 0.5*v4;
          c.d_vatom(i,5) += 0.5*v5;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    compute_item<0,0>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type &energy_virial) const {
    energy_virial += compute_item<1,0>(i,list);
  }

  KOKKOS_INLINE_FUNCTION
  static void init(volatile value_type &update) {
    update.evdwl = 0;
    update.ecoul = 0;
    update.v[0] = 0;
    update.v[1] = 0;
    update.v[2] = 0;
    update.v[3] = 0;
    update.v[4] = 0;
    update.v[5] = 0;
  }
  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &source) {
    update.evdwl += source.evdwl;
    update.ecoul += source.ecoul;
    update.v[0] += source.v[0];
    update.v[1] += source.v[1];
    update.v[2] += source.v[2];
    update.v[3] += source.v[3];
    update.v[4] += source.v[4];
    update.v[5] += source.v[5];
  }


};

template<class PairStyle, class Specialisation>
EV_FLOAT pair_compute (PairStyle* fpair, NeighListKokkos<typename PairStyle::device_type>* list) {
  EV_FLOAT ev;
  if(fpair->atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (fpair->neighflag == FULL) {
      PairComputeFunctor<PairStyle,FULL,false,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == HALFTHREAD) {
      PairComputeFunctor<PairStyle,HALFTHREAD,false,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == HALF) {
      PairComputeFunctor<PairStyle,HALF,false,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == N2) {
      PairComputeFunctor<PairStyle,N2,false,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(fpair->nlocal,ff,ev);
      else Kokkos::parallel_for(fpair->nlocal,ff);
    } else if (fpair->neighflag == FULLCLUSTER) {
      typedef PairComputeFunctor<PairStyle,FULLCLUSTER,false,Specialisation >
        f_type;
      f_type ff(fpair, list);
      #ifdef KOKKOS_HAVE_CUDA
        const int teamsize = Kokkos::Impl::is_same<typename f_type::device_type, Kokkos::Cuda>::value ? 256 : 1;
      #else
        const int teamsize = 1;
      #endif
      const int nteams = (list->inum*f_type::vectorization::increment+teamsize-1)/teamsize;
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(Kokkos::ParallelWorkRequest(nteams,teamsize),ff,ev);
      else Kokkos::parallel_for(Kokkos::ParallelWorkRequest(nteams,teamsize),ff);
    }
  } else {
    if (fpair->neighflag == FULL) {
      PairComputeFunctor<PairStyle,FULL,true,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == HALFTHREAD) {
      PairComputeFunctor<PairStyle,HALFTHREAD,true,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == HALF) {
      PairComputeFunctor<PairStyle,HALF,true,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(list->inum,ff,ev);
      else Kokkos::parallel_for(list->inum,ff);
    } else if (fpair->neighflag == N2) {
      PairComputeFunctor<PairStyle,N2,true,Specialisation >
        ff(fpair, list);
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(fpair->nlocal,ff,ev);
      else Kokkos::parallel_for(fpair->nlocal,ff);
    } else if (fpair->neighflag == FULLCLUSTER) {
      typedef PairComputeFunctor<PairStyle,FULLCLUSTER,true,Specialisation >
        f_type;
      f_type ff(fpair, list);
      #ifdef KOKKOS_HAVE_CUDA
        const int teamsize = Kokkos::Impl::is_same<typename f_type::device_type, Kokkos::Cuda>::value ? 256 : 1;
      #else
        const int teamsize = 1;
      #endif
      const int nteams = (list->inum*f_type::vectorization::increment+teamsize-1)/teamsize;
      if (fpair->eflag || fpair->vflag) Kokkos::parallel_reduce(Kokkos::ParallelWorkRequest(nteams,teamsize),ff,ev);
      else Kokkos::parallel_for(Kokkos::ParallelWorkRequest(nteams,teamsize),ff);
    }
  }
  return ev;
}

}

#endif
#endif

/* ERROR/WARNING messages:

*/
