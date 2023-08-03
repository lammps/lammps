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
   Contributing author: Rodrigo Canales (RWTH Aachen University)
------------------------------------------------------------------------- */

#include "pair_buck_coul_cut_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#define C_FORCE_T typename ForceConst<flt_t>::c_force_t
#define C_ENERGY_T typename ForceConst<flt_t>::c_energy_t
#define C_CUT_T typename ForceConst<flt_t>::c_cut_t

PairBuckCoulCutIntel::PairBuckCoulCutIntel(LAMMPS *lmp) :
  PairBuckCoulCut(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

void PairBuckCoulCutIntel::compute(int eflag, int vflag)
{
  if (fix->precision()==FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  else if (fix->precision()==FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);

  fix->balance_stamp();
  vflag_fdotr = 0;
}

template <class flt_t, class acc_t>
void PairBuckCoulCutIntel::compute(int eflag, int vflag,
                                   IntelBuffers<flt_t,acc_t> *buffers,
                                   const ForceConst<flt_t> &fc)
{
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"INTEL package does not support per-atom stress");
  if (vflag && !vflag_fdotr && force->newton_pair)
    error->all(FLERR,"INTEL package does not support pair_modify nofdotr "
               "with newton on");

  const int inum = list->inum;
  const int nthreads = comm->nthreads;
  const int host_start = fix->host_start_pair();
  const int offload_end = fix->offload_end_pair();
  const int ago = neighbor->ago;

  if (ago != 0 && fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);

    int packthreads;
    if (nthreads > INTEL_HTHREADS) packthreads = nthreads;
    else packthreads = 1;
    #if defined(_OPENMP)
    #pragma omp parallel if (packthreads > 1)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost,
                                packthreads, sizeof(ATOM_T));
      buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  int ovflag = 0;
  if (vflag_fdotr) ovflag = 2;
  else if (vflag) ovflag = 1;
  if (eflag) {
    if (force->newton_pair) {
      eval<1,1>(1, ovflag, buffers, fc, 0, offload_end);
      eval<1,1>(0, ovflag, buffers, fc, host_start, inum);
    } else {
      eval<1,0>(1, ovflag, buffers, fc, 0, offload_end);
      eval<1,0>(0, ovflag, buffers, fc, host_start, inum);
    }
  } else {
    if (force->newton_pair) {
      eval<0,1>(1, ovflag, buffers, fc, 0, offload_end);
      eval<0,1>(0, ovflag, buffers, fc, host_start, inum);
    } else {
      eval<0,0>(1, ovflag, buffers, fc, 0, offload_end);
      eval<0,0>(0, ovflag, buffers, fc, host_start, inum);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairBuckCoulCutIntel::eval(const int offload, const int vflag,
                                IntelBuffers<flt_t,acc_t> *buffers,
                                const ForceConst<flt_t> &fc,
                                const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  flt_t * _noalias const q = buffers->get_q(offload);

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;  // NOLINT

  const flt_t * _noalias const special_coul = fc.special_coul;
  const flt_t * _noalias const special_lj = fc.special_lj;
  const flt_t qqrd2e = force->qqrd2e;

  const C_FORCE_T * _noalias const c_force = fc.c_force[0];
  const C_ENERGY_T * _noalias const c_energy = fc.c_energy[0];
  const C_CUT_T * _noalias const c_cut = fc.c_cut[0];

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, NEWTON_PAIR, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);

  const int nthreads = tc;
  #ifdef _LMP_INTEL_OFFLOAD
  int *overflow = fix->get_off_overflow_flag();
  double *timer_compute = fix->off_watch_pair();
  // Redeclare as local variables for offload
  const int ncoulmask = this->ncoulmask;
  const int ncoulshiftbits = this->ncoulshiftbits;

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if (offload)                 \
    in(special_lj,special_coul:length(0) alloc_if(0) free_if(0)) \
    in(c_force, c_energy, c_cut:length(0) alloc_if(0) free_if(0))      \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(q:length(q_size) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(astart,nthreads,qqrd2e,inum,nall,ntypes,vflag,eatom) \
    in(f_stride,nlocal,minlocal,separate_flag,offload) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(NEWTON_PAIR, separate_flag, nlocal, nall,
                              f_stride, x, q);

    acc_t oevdwl, oecoul, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = oecoul = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    if (NEWTON_PAIR == 0 && inum != nlocal)
      memset(f_start, 0, f_stride * sizeof(FORCE_T));

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,oecoul,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iip, iito, tid;
      IP_PRE_omp_stride_id(iifrom, iip, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      int foff;
      if (NEWTON_PAIR) foff = tid * f_stride - minlocal;
      else foff = -minlocal;
      FORCE_T * _noalias const f = f_start + foff;
      if (NEWTON_PAIR) memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        const int itype = x[i].w;

        const int ptr_off = itype * ntypes;
        const C_FORCE_T * _noalias const c_forcei = c_force + ptr_off;
        const C_ENERGY_T * _noalias const c_energyi = c_energy + ptr_off;
        const C_CUT_T * _noalias const c_cuti = c_cut + ptr_off;
        const int   * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum, offload);

        acc_t fxtmp,fytmp,fztmp,fwtmp;
        acc_t sevdwl, secoul, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const flt_t qtmp = q[i];
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EFLAG) fwtmp = sevdwl = secoul = (acc_t)0;
        if (NEWTON_PAIR == 0)
          if (vflag == VIRIAL_PAIR)
            sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                                   sv0, sv1, sv2, sv3, sv4, sv5)
#else
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                               sv0, sv1, sv2, sv3, sv4, sv5)
#endif
        #pragma vector aligned
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          flt_t forcecoul, forcebuck, evdwl, ecoul;
          forcecoul = forcebuck = evdwl = ecoul = (flt_t)0.0;

          const int sbindex = jlist[jj] >> SBBITS & 3;
          const int j = jlist[jj] & NEIGHMASK;

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const int jtype = IP_PRE_dword_index(x[j].w);
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          const flt_t r = std::sqrt(rsq);
          const flt_t r2inv = (flt_t)1.0 / rsq;

          #ifdef INTEL_VMASK
          if (rsq < c_cuti[jtype].cut_coulsq) {
          #endif
            forcecoul = qqrd2e * qtmp*q[j]/r;
            if (EFLAG)
              ecoul = forcecoul;
            if (sbindex) {
              const flt_t factor_coul = special_coul[sbindex];
              forcecoul *= factor_coul;
              if (EFLAG)
                ecoul *= factor_coul;

            }
          #ifdef INTEL_VMASK
          }
          #else
          if (rsq >= c_cuti[jtype].cut_coulsq)
            { forcecoul = (flt_t)0.0; ecoul = (flt_t)0.0; }
          #endif

          #ifdef INTEL_VMASK
          if (rsq < c_cuti[jtype].cut_ljsq) {
          #endif
            flt_t r6inv = r2inv * r2inv * r2inv;
            flt_t rexp = std::exp(-r * c_forcei[jtype].rhoinv);
            forcebuck = r * rexp * c_forcei[jtype].buck1 -
              r6inv * c_forcei[jtype].buck2;
            if (EFLAG)
              evdwl = rexp * c_energyi[jtype].a -
                r6inv * c_energyi[jtype].c -
                c_energyi[jtype].offset;
            if (sbindex) {
              const flt_t factor_lj = special_lj[sbindex];
              forcebuck *= factor_lj;
              if (EFLAG)
                evdwl *= factor_lj;
            }
          #ifdef INTEL_VMASK
          }
          #else
          if (rsq >= c_cuti[jtype].cut_ljsq)
            { forcebuck = (flt_t)0.0; evdwl = (flt_t)0.0; }
          #endif

          #ifdef INTEL_VMASK
          if (rsq < c_cuti[jtype].cutsq) {
          #endif
            const flt_t fpair = (forcecoul + forcebuck) * r2inv;
            const flt_t fpx = fpair * delx;
            fxtmp += fpx;
            if (NEWTON_PAIR) f[j].x -= fpx;
            const flt_t fpy = fpair * dely;
            fytmp += fpy;
            if (NEWTON_PAIR) f[j].y -= fpy;
            const flt_t fpz = fpair * delz;
            fztmp += fpz;
            if (NEWTON_PAIR) f[j].z -= fpz;


            if (EFLAG) {
              sevdwl += evdwl;
              secoul += ecoul;
              if (eatom) {
                fwtmp += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
                if (NEWTON_PAIR)
                  f[j].w += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
              }
            }
            if (NEWTON_PAIR == 0)
              IP_PRE_ev_tally_nborv(vflag, delx, dely, delz, fpx, fpy, fpz);
          #ifdef INTEL_VMASK
          }
          #endif
        } // for jj
        if (NEWTON_PAIR) {
          f[i].x += fxtmp;
          f[i].y += fytmp;
          f[i].z += fztmp;
        } else {
          f[i].x = fxtmp;
          f[i].y = fytmp;
          f[i].z = fztmp;
        }
        IP_PRE_ev_tally_atomq(NEWTON_PAIR, EFLAG, vflag, f, fwtmp);
      } // for ii

      IP_PRE_fdotr_reduce_omp(NEWTON_PAIR, nall, minlocal, nthreads, f_start,
                              f_stride, x, offload, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
    } // end of omp parallel region

    IP_PRE_fdotr_reduce(NEWTON_PAIR, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      if (NEWTON_PAIR == 0) {
        oevdwl *= (acc_t)0.5;
        oecoul *= (acc_t)0.5;
        ov0 *= (acc_t)0.5;
        ov1 *= (acc_t)0.5;
        ov2 *= (acc_t)0.5;
        ov3 *= (acc_t)0.5;
        ov4 *= (acc_t)0.5;
        ov5 *= (acc_t)0.5;
      }
      ev_global[0] = oevdwl;
      ev_global[1] = oecoul;
      ev_global[2] = ov0;
      ev_global[3] = ov1;
      ev_global[4] = ov2;
      ev_global[5] = ov3;
      ev_global[6] = ov4;
      ev_global[7] = ov5;
    }
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end of offload region

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, nullptr, offload);
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulCutIntel::init_style()
{
  PairBuckCoulCut::init_style();
  if (force->newton_pair == 0)
    neighbor->find_request(this)->enable_full();

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = fix->coprocessor_number();
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());

}

template <class flt_t, class acc_t>
void PairBuckCoulCutIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;
  int ntable = 1;
  if (ncoultablebits)
    for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  fc.set_ntypes(tp1, ntable, memory, _cop);

  // Repeat cutsq calculation because done after call to init_style
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      double cut;
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0))
        cut = init_one(i, j);
      else
        cut = 0.0;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_coul[i] = force->special_coul[i];
    fc.special_coul[0] = 1.0;
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.c_cut[i][j].cutsq = cutsq[i][j];
      fc.c_cut[i][j].cut_ljsq = cut_ljsq[i][j];
      fc.c_cut[i][j].cut_coulsq = cut_coulsq[i][j];
      fc.c_force[i][j].buck1 = buck1[i][j];
      fc.c_force[i][j].buck2 = buck2[i][j];
      fc.c_force[i][j].rhoinv = rhoinv[i][j];
      fc.c_energy[i][j].a = a[i][j];
      fc.c_energy[i][j].c = c[i][j];
      fc.c_energy[i][j].offset = offset[i][j];
    }
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  flt_t * special_lj = fc.special_lj;
  flt_t * special_coul = fc.special_coul;
  C_FORCE_T * c_force = fc.c_force[0];
  C_ENERGY_T * c_energy = fc.c_energy[0];
  C_CUT_T * c_cut = fc.c_cut[0];
  int tp1sq = tp1 * tp1;
  #pragma offload_transfer target(mic:_cop) \
    in(special_lj, special_coul: length(4) alloc_if(0) free_if(0)) \
    in(c_force, c_energy, c_cut: length(tp1sq) alloc_if(0) free_if(0))
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairBuckCoulCutIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                           const int ntable,
                                                           Memory *memory,
                                                           const int cop) {
  if (memory != nullptr) _memory = memory;
  if ((ntypes != _ntypes) || (ntable != _ntable)) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      c_force_t * oc_force = c_force[0];
      c_energy_t * oc_energy = c_energy[0];
      c_cut_t * oc_cut = c_cut[0];

      if (ospecial_lj != nullptr && oc_force != nullptr && oc_cut != nullptr &&
          oc_energy != nullptr && ospecial_coul != nullptr &&
          _cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj, ospecial_coul: alloc_if(0) free_if(1)) \
          nocopy(oc_force, oc_energy: alloc_if(0) free_if(1))        \
          nocopy(oc_cut: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(c_force);
      _memory->destroy(c_energy);
      _memory->destroy(c_cut);
    }
    if (ntypes > 0) {
      _cop = cop;
      _memory->create(c_force,ntypes,ntypes,"fc.c_force");
      _memory->create(c_energy,ntypes,ntypes,"fc.c_energy");
      _memory->create(c_cut,ntypes,ntypes,"fc.c_cut");


      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      c_force_t * oc_force = c_force[0];
      c_energy_t * oc_energy = c_energy[0];
      c_cut_t * oc_cut = c_cut[0];
      int tp1sq = ntypes*ntypes;
      if (ospecial_lj != nullptr && oc_force != nullptr && oc_cut != nullptr &&
          oc_energy != nullptr && ospecial_coul != nullptr &&
          cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj: length(4) alloc_if(1) free_if(0)) \
          nocopy(ospecial_coul: length(4) alloc_if(1) free_if(0)) \
          nocopy(oc_force: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_energy: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_cut: length(tp1sq) alloc_if(1) free_if(0))

      }
      #endif
    }
  }
  _ntypes=ntypes;
  _ntable=ntable;
}
