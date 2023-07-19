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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "pair_sw_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"

#include <cstring>

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif
#include <cmath>
#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#ifdef LMP_USE_AVXCD
#include "intel_simd.h"
using namespace ip_simd;
#endif

using namespace LAMMPS_NS;

#define FC_PACKED0_T typename ForceConst<flt_t>::fc_packed0
#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED1p2_T typename ForceConst<flt_t>::fc_packed1p2
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2
#define FC_PACKED3_T typename ForceConst<flt_t>::fc_packed3

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWIntel::PairSWIntel(LAMMPS *lmp) : PairSW(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void PairSWIntel::compute(int eflag, int vflag)
{
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);

  fix->balance_stamp();
  vflag_fdotr = 0;
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairSWIntel::compute(int eflag, int vflag,
                          IntelBuffers<flt_t,acc_t> *buffers,
                          const ForceConst<flt_t> &fc)
{
  ev_init(eflag, vflag);
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
      buffers->thr_pack(ifrom, ito, ago);
    }

    fix->stop_watch(TIME_PACK);
  }

  int ovflag = 0;
  if (vflag_fdotr) ovflag = 2;
  else if (vflag) ovflag = 1;
  if (_onetype) {
    if (_spq) {
      if (eflag) {
        eval<1,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (eflag) {
        eval<0,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  } else {
    if (_spq) {
      if (eflag) {
        eval<1,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (eflag) {
        eval<0,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
#ifndef LMP_USE_AVXCD

template <int SPQ,int ONETYPE,int EFLAG,class flt_t,class acc_t>
void PairSWIntel::eval(const int offload, const int vflag,
                       IntelBuffers<flt_t,acc_t> *buffers,
                       const ForceConst<flt_t> &fc, const int astart,
                       const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * _noalias const firstneigh = list->firstneigh[ilist[0]];

  int *nhalf, *cnum;
  buffers->get_list_data3(list, nhalf, cnum);
  const int * _noalias const numneighhalf = nhalf;
  const int * _noalias const cnumneigh = cnum;

  const FC_PACKED0_T * _noalias const p2 = fc.p2[0];
  const FC_PACKED1_T * _noalias const p2f = fc.p2f[0];
  const FC_PACKED1p2_T * _noalias const p2f2 = fc.p2f2[0];
  const FC_PACKED2_T * _noalias const p2e = fc.p2e[0];
  const FC_PACKED3_T * _noalias const p3 = fc.p3[0][0];

  flt_t * _noalias const ccachex = buffers->get_ccachex();
  flt_t * _noalias const ccachey = buffers->get_ccachey();
  flt_t * _noalias const ccachez = buffers->get_ccachez();
  flt_t * _noalias const ccachew = buffers->get_ccachew();
  int * _noalias const ccachei = buffers->get_ccachei();
  int * _noalias const ccachej = buffers->get_ccachej();
  const int ccache_stride = _ccache_stride;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;
  const int onetype = _onetype;
  const int onetype3 = _onetype3;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, /* NEWTON_PAIR*/ 1, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  #ifdef _LMP_INTEL_OFFLOAD
  double *timer_compute = fix->off_watch_pair();
  int *overflow = fix->get_off_overflow_flag();

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if (offload) \
    in(p2,p2f,p2f2,p2e,p3:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(ccachex,ccachey,ccachez,ccachew:length(0) alloc_if(0) free_if(0)) \
    in(ccachei,ccachej:length(0) alloc_if(0) free_if(0)) \
    in(ccache_stride,nthreads,inum,nall,ntypes,vflag,eatom,offload,onetype3) \
    in(astart,nlocal,f_stride,minlocal,separate_flag,onetype) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(1, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iip, iito, tid;
      IP_PRE_omp_stride_id(iifrom, iip, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;

      // loop over full neighbor list of my atoms
      flt_t cutsq, cut, powerp, powerq, sigma, c1, c2, c3, c4, c5, c6;
      flt_t sigma_gamma, costheta, lambda_epsilon, lambda_epsilon2;
      if (ONETYPE) {
        cutsq = p2[onetype].cutsq;
        cut = p2f[onetype].cut;
        sigma = p2f[onetype].sigma;
        c1 = p2f2[onetype].c1;
        c2 = p2f2[onetype].c2;
        c3 = p2f2[onetype].c3;
        c4 = p2f2[onetype].c4;
        sigma_gamma = p2[onetype].sigma_gamma;
        costheta = p3[onetype3].costheta;
        lambda_epsilon = p3[onetype3].lambda_epsilon;
        lambda_epsilon2 = p3[onetype3].lambda_epsilon2;
        if (SPQ == 0) {
          powerp = p2f[onetype].powerp;
          powerq = p2f[onetype].powerq;
        }
        if (EFLAG) {
          c5 = p2e[onetype].c5;
          c6 = p2e[onetype].c6;
        }
      }

      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        int itype, itype_offset;
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;

        if (!ONETYPE) {
          itype = x[i].w;
          itype_offset = itype * ntypes;
        }

        const int * _noalias const jlist = firstneigh + cnumneigh[ii];
        const int jnum = numneigh[i];
        const int jnumhalf = numneighhalf[ii];

        acc_t fxtmp, fytmp, fztmp, fwtmp;
        acc_t sevdwl;
        fxtmp = fytmp = fztmp = (acc_t)0.0;
        if (EFLAG) fwtmp = sevdwl = (acc_t)0;

        int ejnum = 0, ejnumhalf = 0;
        #pragma vector aligned
        #pragma ivdep
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;
          const flt_t delx = x[j].x - xtmp;
          const flt_t dely = x[j].y - ytmp;
          const flt_t delz = x[j].z - ztmp;
          int jtype, ijtype;
          if (!ONETYPE) {
            jtype = x[j].w;
            ijtype = IP_PRE_dword_index(itype_offset + jtype);
            cutsq = p2[ijtype].cutsq;
          }
          const flt_t rsq1 = delx * delx + dely * dely + delz * delz;
          if (rsq1 < cutsq) {
            tdelx[ejnum] = delx;
            tdely[ejnum] = dely;
            tdelz[ejnum] = delz;
            trsq[ejnum] = rsq1;
            tj[ejnum] = j;
            if (!ONETYPE) tjtype[ejnum] = jtype;
            ejnum++;
            if (jj < jnumhalf) ejnumhalf++;
          }
        }

        int ejnum_pad = ejnum;
        IP_PRE_neighbor_pad(ejnum_pad, offload);
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma loop_count min=1, max=INTEL_COMPILE_WIDTH-1, \
                avg=INTEL_COMPILE_WIDTH/2
        #endif
        for (int jj = ejnum; jj < ejnum_pad; jj++) {
          tdelx[jj] = (flt_t)0.0;
          tdely[jj] = (flt_t)0.0;
          tdelz[jj] = (flt_t)0.0;
          trsq[jj] = p2[3].cutsq + (flt_t)1.0;
          tj[jj] = nall;
          if (!ONETYPE) tjtype[jj] = 0;
        }

        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl)
#else
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl)
#endif
        #pragma vector aligned
        #endif
        for (int jj = 0; jj < ejnum_pad; jj++) {
          acc_t fjxtmp, fjytmp, fjztmp, fjtmp;
          fjxtmp = fjytmp = fjztmp = (acc_t)0.0;
          if (EFLAG) fjtmp = (acc_t)0.0;
          int ijtype;

          if (!ONETYPE) ijtype = IP_PRE_dword_index(tjtype[jj] + itype_offset);
          const flt_t rsq1 = trsq[jj];

          const flt_t rinvsq1 = (flt_t)1.0 / rsq1;
          const flt_t r1 = (flt_t)1.0/std::sqrt(rinvsq1);
          if (!ONETYPE) cut = p2f[ijtype].cut;
          const flt_t rainv1 = (flt_t)1.0 / (r1 - cut);

          // two-body interactions, skip half of them
          flt_t rp, rq;
          if (SPQ == 1) {
            rp = r1 * r1;
            rp *= rp;
            rp = (flt_t)1.0 / rp;
            rq = (flt_t)1.0;
          } else {
            if (!ONETYPE) {
              powerp = p2f[ijtype].powerp;
              powerq = p2f[ijtype].powerq;
            }
            rp = std::pow(r1, powerp);
            rq = std::pow(r1, powerq);
          }

          if (!ONETYPE) {
            sigma = p2f[ijtype].sigma;
            c1 = p2f2[ijtype].c1;
            c2 = p2f2[ijtype].c2;
            c3 = p2f2[ijtype].c3;
            c4 = p2f2[ijtype].c4;
          }

          const flt_t rainvsq = rainv1 * rainv1 * r1;
          flt_t expsrainv = std::exp(sigma * rainv1);
          if (jj >= ejnumhalf) expsrainv = (flt_t)0.0;
          const flt_t fpair = (c1 * rp - c2 * rq + (c3 * rp - c4 * rq) *
                               rainvsq) * expsrainv * rinvsq1;

          const flt_t delx = tdelx[jj];
          const flt_t dely = tdely[jj];
          const flt_t delz = tdelz[jj];
          const flt_t fpx = fpair * delx;
          fxtmp -= fpx;
          fjxtmp += fpx;
          const flt_t fpy = fpair * dely;
          fytmp -= fpy;
          fjytmp += fpy;
          const flt_t fpz = fpair * delz;
          fztmp -= fpz;
          fjztmp += fpz;

          if (EFLAG) {
            flt_t evdwl;
            if (!ONETYPE) {
              c5 = p2e[ijtype].c5;
              c6 = p2e[ijtype].c6;
            }
            evdwl = (c5 * rp - c6 * rq) * expsrainv;
            sevdwl += evdwl;
            if (eatom) {
              fwtmp += (flt_t)0.5 * evdwl;
              fjtmp += (flt_t)0.5 * evdwl;
            }
          }

          /*---------------------------------------------*/

          int ijkoff;
          if (!ONETYPE) {
            sigma_gamma = p2[ijtype].sigma_gamma;
            ijkoff = ijtype * ntypes;
          }

          flt_t gsrainv1 = sigma_gamma * rainv1;
          flt_t gsrainvsq1 = gsrainv1 * rainv1 / r1;
          flt_t expgsrainv1 = std::exp(gsrainv1);

          for (int kk = 0; kk < ejnum; kk++) {
            int iktype, ijktype;
            if (!ONETYPE) {
              iktype = tjtype[kk];
              ijktype = IP_PRE_dword_index(ijkoff + iktype);
              iktype = IP_PRE_dword_index(iktype + itype_offset);
              cut = p2[iktype].cut;
              sigma_gamma = p2[iktype].sigma_gamma;
              costheta = p3[ijktype].costheta;
              lambda_epsilon = p3[ijktype].lambda_epsilon;
              lambda_epsilon2 = p3[ijktype].lambda_epsilon2;
            }

            flt_t delr2[3];
            delr2[0] = tdelx[kk];
            delr2[1] = tdely[kk];
            delr2[2] = tdelz[kk];
            const flt_t rsq2 = trsq[kk];

            const flt_t rinvsq2 = (flt_t)1.0 / rsq2;
            const flt_t r2 = (flt_t)1.0 / std::sqrt(rinvsq2);
            const flt_t rainv2 = (flt_t)1.0 / (r2 - cut);
            const flt_t gsrainv2 = sigma_gamma * rainv2;
            const flt_t gsrainvsq2 = gsrainv2 * rainv2 / r2;
            const flt_t expgsrainv2 = std::exp(gsrainv2);

            const flt_t rinv12 = (flt_t)1.0 / (r1 * r2);
            const flt_t cs = (delx * delr2[0] + dely * delr2[1] +
                              delz * delr2[2]) * rinv12;
            const flt_t delcs = cs - costheta;
            const flt_t delcssq = delcs*delcs;

            flt_t kfactor;
            if (jj == kk || jj >= ejnum) kfactor = (flt_t)0.0;
            else kfactor = (flt_t)1.0;

            const flt_t facexp = expgsrainv1*expgsrainv2*kfactor;
            const flt_t facrad = lambda_epsilon * facexp * delcssq;
            const flt_t frad1 = facrad*gsrainvsq1;
            const flt_t frad2 = facrad*gsrainvsq2;
            const flt_t facang = lambda_epsilon2 * facexp * delcs;
            const flt_t facang12 = rinv12*facang;
            const flt_t csfacang = cs*facang;
            const flt_t csfac1 = rinvsq1*csfacang;

            const flt_t fjx = delx*(frad1+csfac1)-delr2[0]*facang12;
            const flt_t fjy = dely*(frad1+csfac1)-delr2[1]*facang12;
            const flt_t fjz = delz*(frad1+csfac1)-delr2[2]*facang12;

            fxtmp -= fjx;
            fytmp -= fjy;
            fztmp -= fjz;
            fjxtmp += fjx;
            fjytmp += fjy;
            fjztmp += fjz;

            if (EFLAG) {
              const flt_t evdwl = facrad * (flt_t)0.5;
              sevdwl += evdwl;
              if (eatom) {
                fwtmp += (acc_t)0.33333333 * evdwl;
                fjtmp += (acc_t)0.33333333 * facrad;
              }
            }
          } // for kk
          const int j = IP_PRE_dword_index(tj[jj]);
          f[j].x += fjxtmp;
          f[j].y += fjytmp;
          f[j].z += fjztmp;
          if (EFLAG)
            if (eatom) f[j].w += fjtmp;
        } // for jj

        f[i].x += fxtmp;
        f[i].y += fytmp;
        f[i].z += fztmp;

        if (EFLAG) {
          f[i].w += fwtmp;
          oevdwl += sevdwl;
        }
      } // for ii

      IP_PRE_fdotr_reduce_omp(1, nall, minlocal, nthreads, f_start, f_stride,
                              x, offload, vflag, ov0, ov1, ov2, ov3, ov4, ov5);
    } // end omp

    IP_PRE_fdotr_reduce(1, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      ev_global[0] = oevdwl;
      ev_global[1] = (acc_t)0.0;
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
  } // end offload
  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, nullptr, offload);
}

#else

/* ----------------------------------------------------------------------

Vector intrinsics are temporarily being used for the Stillinger-Weber
potential to allow for advanced features in the AVX512 instruction set to
be exploited on early hardware. We hope to see compiler improvements for
AVX512 that will eliminate this requirement, so it is not recommended to
develop code based on the intrinsics implementation. Please e-mail the
authors for more details.

------------------------------------------------------------------------- */

template <int SPQ,int ONETYPE,int EFLAG,class flt_t,class acc_t>
void PairSWIntel::eval(const int offload, const int vflag,
                       IntelBuffers<flt_t,acc_t> *buffers,
                       const ForceConst<flt_t> &fc, const int astart,
                       const int aend)
{
  typedef typename SIMD_type<flt_t>::SIMD_vec SIMD_flt_t;
  typedef typename SIMD_type<acc_t>::SIMD_vec SIMD_acc_t;
  const int swidth = SIMD_type<flt_t>::width();

  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * _noalias const firstneigh = list->firstneigh[ilist[0]];

  int *nhalf, *cnum;
  buffers->get_list_data3(list, nhalf, cnum);
  const int * _noalias const numneighhalf = nhalf;
  const int * _noalias const cnumneigh = cnum;

  const FC_PACKED0_T * _noalias const p2 = fc.p2[0];
  const FC_PACKED1_T * _noalias const p2f = fc.p2f[0];
  const FC_PACKED1p2_T * _noalias const p2f2 = fc.p2f2[0];
  const FC_PACKED2_T * _noalias const p2e = fc.p2e[0];
  const FC_PACKED3_T * _noalias const p3 = fc.p3[0][0];

  flt_t * _noalias const ccachex = buffers->get_ccachex();
  flt_t * _noalias const ccachey = buffers->get_ccachey();
  flt_t * _noalias const ccachez = buffers->get_ccachez();
  flt_t * _noalias const ccachew = buffers->get_ccachew();
  int * _noalias const ccachei = buffers->get_ccachei();
  int * _noalias const ccachej = buffers->get_ccachej();
  acc_t * _noalias const ccachef = buffers->get_ccachef();
  const int ccache_stride = _ccache_stride;
  const int ccache_stride3 = _ccache_stride3;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, /* NEWTON_PAIR*/ 1, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  #ifdef _LMP_INTEL_OFFLOAD
  double *timer_compute = fix->off_watch_pair();
  int *overflow = fix->get_off_overflow_flag();

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if (offload) \
    in(p2,p2f,p2f2,p2e,p3:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(ccachex,ccachey,ccachez,ccachew:length(0) alloc_if(0) free_if(0)) \
    in(ccachei,ccachej,ccachef:length(0) alloc_if(0) free_if(0)) \
    in(ccache_stride,nthreads,inum,nall,ntypes,vflag,eatom,offload) \
    in(astart,nlocal,f_stride,minlocal,separate_flag) \
    in(ccache_stride3)                                          \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(1, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iip, iito, tid;
      IP_PRE_omp_stride_id_vec(iifrom, iip, iito, tid, inum, nthreads,
                               swidth);

      iifrom += astart;
      iito += astart;

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;
      acc_t * _noalias const tf = ccachef + tid * ccache_stride3;

      // loop over full neighbor list of my atoms

      SIMD_flt_t cutsq, cut, powerp, powerq, sigma, c1, c2, c3,c4, c5, c6;
      SIMD_flt_t sigma_gamma, costheta, lambda_epsilon, lambda_epsilon2;
      if (ONETYPE) {
        cutsq = SIMD_set(p2[_onetype].cutsq);
        cut = SIMD_set(p2f[_onetype].cut);
        sigma = SIMD_set(p2f[_onetype].sigma);
        c1 = SIMD_set(p2f2[_onetype].c1);
        c2 = SIMD_set(p2f2[_onetype].c2);
        c3 = SIMD_set(p2f2[_onetype].c3);
        c4 = SIMD_set(p2f2[_onetype].c4);
        sigma_gamma = SIMD_set(p2[_onetype].sigma_gamma);
        costheta = SIMD_set(p3[_onetype3].costheta);
        lambda_epsilon = SIMD_set(p3[_onetype3].lambda_epsilon);
        lambda_epsilon2 = SIMD_set(p3[_onetype3].lambda_epsilon2);
        if (SPQ == 0) {
          powerp = SIMD_set(p2f[_onetype].powerp);
          powerq = SIMD_set(p2f[_onetype].powerq);
        }
        if (EFLAG) {
          c5 = SIMD_set(p2e[_onetype].c5);
          c6 = SIMD_set(p2e[_onetype].c6);
        }
      }

      acc_t * const dforce = &(f[0].x);
      for (int i = iifrom; i < iito; i += iip) {
        SIMD_int ilistv = SIMD_load(ilist + i);
        SIMD_int goffset = ilistv * 16;
        SIMD_mask imask;
        if (swidth == 16)
          imask = 0xFFFF;
        else
          imask = 0xFF;
        const int rem = iito - i;
        if (rem < swidth) imask = imask >> (swidth - rem);
        SIMD_flt_t xtmp, ytmp, ztmp;
        SIMD_int itype, itype_offset;

        if (ONETYPE)
          SIMD_atom_gather(imask, &(x[0].x), goffset, xtmp, ytmp, ztmp);
        else {
          SIMD_atom_gather(imask, &(x[0].x), goffset, xtmp, ytmp, ztmp, itype);
          itype_offset = itype * ntypes;
        }

        #ifdef LMP_INTEL_3BODY_FAST
        const int* ng = firstneigh + cnumneigh[i] - swidth;
        const SIMD_int jnum = SIMD_loadz(imask, numneigh + i);
        #else
        SIMD_int ng = SIMD_load(cnumneigh + i);
        const SIMD_int jnum = SIMD_gatherz(imask, numneigh, ilistv);
        ng = ng - 1;
        #endif
        const SIMD_int jnumhalf = SIMD_loadz(imask, numneighhalf + i);
        const int jnum_max = SIMD_max(jnum);

        SIMD_acc_t fxtmp = SIMD_set((acc_t)0);
        SIMD_acc_t fytmp = SIMD_set((acc_t)0);
        SIMD_acc_t fztmp = SIMD_set((acc_t)0);
        SIMD_acc_t fwtmp, fxtmp2, fytmp2, fztmp2, fwtmp2;
        if (is_same<flt_t,acc_t>::value == 0) {
          fxtmp2 = SIMD_set((acc_t)0);
          fytmp2 = SIMD_set((acc_t)0);
          fztmp2 = SIMD_set((acc_t)0);
          if (EFLAG) fwtmp2 = SIMD_set((acc_t)0);
        }

        SIMD_acc_t sevdwl;
        if (EFLAG) {
          fwtmp = SIMD_set((acc_t)0);
          sevdwl = SIMD_set((acc_t)0);
        }

        SIMD_int ejnum = SIMD_set(0);
        SIMD_int ejnumhalf = SIMD_set(0);
        SIMD_int coffset = SIMD_set(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                    11, 12, 13, 14, 15);
        for (int jj = 0; jj < jnum_max; jj++) {
          SIMD_mask jmask = jj < jnum;

          #ifdef LMP_INTEL_3BODY_FAST
          ng += swidth;
          SIMD_int j = SIMD_load(ng);
          #else
          ng = ng + 1;
          SIMD_int j = SIMD_gather(jmask, firstneigh, ng);
          #endif
          j = j & SIMD_set(NEIGHMASK);
          const SIMD_int joffset = j << 4;

          SIMD_flt_t delx, dely, delz;
          SIMD_int jtype, ijtype;
          if (ONETYPE)
            SIMD_atom_gather(jmask, &(x[0].x), joffset, delx, dely, delz);
          else {
            SIMD_atom_gather(jmask, &(x[0].x), joffset, delx, dely, delz,
                             jtype);
            ijtype = (jtype + itype_offset) << 2;
            cutsq = SIMD_gather(jmask, &(p2[0].cutsq), ijtype);
          }

          delx = delx - xtmp;
          dely = dely - ytmp;
          delz = delz - ztmp;
          SIMD_flt_t rsq1 = delx * delx;
          rsq1 = SIMD_fma(dely, dely, rsq1);
          rsq1 = SIMD_fma(delz, delz, rsq1);

          const SIMD_mask rmask = SIMD_lt(jmask, rsq1, cutsq);
          SIMD_scatter(rmask, tdelx, coffset, delx);
          SIMD_scatter(rmask, tdely, coffset, dely);
          SIMD_scatter(rmask, tdelz, coffset, delz);
          SIMD_scatter(rmask, trsq, coffset, rsq1);
          SIMD_scatter(rmask, tj, coffset, j);
          if (!ONETYPE) SIMD_scatter(rmask, tjtype, coffset, jtype);
          ejnum = SIMD_add(rmask, ejnum, 1);
          coffset = SIMD_add(rmask, coffset, swidth);
          const SIMD_mask hmask = SIMD_lt(rmask, SIMD_set(jj), jnumhalf);
          ejnumhalf = SIMD_add(hmask, ejnumhalf, 1);
        }

        const int ejnum_max = SIMD_max(ejnum);
        const int ejnumhalf_max = SIMD_max(ejnumhalf);
        memset(tf, 0, ejnum_max * sizeof(acc_t) * swidth * 3);
        for (int jj = 0; jj < ejnum_max; jj++) {
          SIMD_int ijtype;
          const int coffset = jj * swidth;
          if (!ONETYPE) {
            ijtype = SIMD_load(tjtype + coffset);
            ijtype = (ijtype + itype_offset) << 2;
            cut = SIMD_gather(&(p2f[0].cut), ijtype);
          }

          SIMD_acc_t fjxtmp = SIMD_set((acc_t)0);
          SIMD_acc_t fjytmp = SIMD_set((acc_t)0);
          SIMD_acc_t fjztmp = SIMD_set((acc_t)0);
          SIMD_acc_t fjtmp, fjxtmp2, fjytmp2, fjztmp2, fjtmp2;
          if (EFLAG) fjtmp = SIMD_set((acc_t)0.0);

          if (is_same<flt_t,acc_t>::value == 0) {
            fjxtmp2 = SIMD_set((acc_t)0);
            fjytmp2 = SIMD_set((acc_t)0);
            fjztmp2 = SIMD_set((acc_t)0);
            if (EFLAG) fjtmp2 = SIMD_set((acc_t)0.0);
          }

          const SIMD_flt_t delx = SIMD_load(tdelx + coffset);
          const SIMD_flt_t dely = SIMD_load(tdely + coffset);
          const SIMD_flt_t delz = SIMD_load(tdelz + coffset);
          const SIMD_flt_t rsq1 = SIMD_load(trsq + coffset);

          const SIMD_flt_t rinvsq1 = SIMD_rcp(rsq1);
          const SIMD_flt_t r1 = SIMD_invsqrt(rinvsq1);
          const SIMD_flt_t rainv1 = SIMD_rcp(r1 - cut);

          // two-body interactions, skip half of them
          if (jj < ejnumhalf_max) {
            SIMD_flt_t rp, rq;
            if (SPQ == 1) {
              rp = r1 * r1;
              rp = rp * rp;
              rp = SIMD_rcp(rp);
              rq = SIMD_set((flt_t)1.0);
            } else {
              if (!ONETYPE) {
                powerp = SIMD_gather(&(p2f[0].powerp), ijtype);
                powerq = SIMD_gather(&(p2f[0].powerq), ijtype);
              }
              rp = SIMD_pow(r1, powerp);
              rq = SIMD_pow(r1, powerq);
            }

            if (!ONETYPE) {
              sigma = SIMD_gather(&(p2f[0].sigma), ijtype);
              c1 = SIMD_gather(&(p2f2[0].c1), ijtype);
              c2 = SIMD_gather(&(p2f2[0].c2), ijtype);
              c3 = SIMD_gather(&(p2f2[0].c3), ijtype);
              c4 = SIMD_gather(&(p2f2[0].c4), ijtype);
            }

            const SIMD_flt_t rainvsq = rainv1 * rainv1 * r1;
            const SIMD_flt_t expsrainv = SIMD_exp(sigma * rainv1);
            const SIMD_flt_t fpair = (c1 * rp - c2 * rq + (c3 * rp - c4 * rq) *
                                      rainvsq) * expsrainv * rinvsq1;

            const SIMD_flt_t fjx = delx * fpair;
            const SIMD_flt_t fjy = dely * fpair;
            const SIMD_flt_t fjz = delz * fpair;

            const SIMD_mask hmask = jj < ejnumhalf;
            SIMD_accumulate3(hmask, fjx, fjy, fjz, fxtmp, fytmp, fztmp,
                             fjxtmp, fjytmp, fjztmp, fxtmp2, fytmp2,
                             fztmp2, fjxtmp2, fjytmp2, fjztmp2);

            if (EFLAG) {
              if (!ONETYPE) {
                c5 = SIMD_gather(&(p2e[0].c5), ijtype);
                c6 = SIMD_gather(&(p2e[0].c6), ijtype);
              }
              SIMD_flt_t evdwl;
              evdwl = (c5 * rp - c6 * rq) * expsrainv;
              SIMD_acc_energy3(hmask, evdwl, eatom, sevdwl, fwtmp, fjtmp,
                               fwtmp2, fjtmp2);
            }
          }

          /*---------------------------------------------*/
          SIMD_int ijkoff;
          if (!ONETYPE) {
            sigma_gamma = SIMD_gather(&(p2[0].sigma_gamma), ijtype);
            ijkoff = ijtype * ntypes;
          }
          const SIMD_flt_t gsrainv1 = sigma_gamma * rainv1;
          const SIMD_flt_t gsrainvsq1 = gsrainv1 * rainv1 / r1;
          const SIMD_flt_t expgsrainv1 = SIMD_exp(gsrainv1);

          const SIMD_mask jmask = jj < ejnum;
          for (int kk = jj+1; kk < ejnum_max; kk++) {
            SIMD_int iktype, ijktype;
            const int kcoffset = kk * swidth;
            if (!ONETYPE) {
              iktype = SIMD_load(tjtype + kcoffset);
              ijktype = ijkoff + (iktype << 2);
              iktype = (iktype + itype_offset) << 2;
              cut = SIMD_gather(&(p2[0].cut), iktype);
              sigma_gamma = SIMD_gather(&(p2[0].sigma_gamma), iktype);
              costheta = SIMD_gather(&(p3[0].costheta), ijktype);
              lambda_epsilon = SIMD_gather(&(p3[0].lambda_epsilon), ijktype);
              lambda_epsilon2 = SIMD_gather(&(p3[0].lambda_epsilon2), ijktype);
            }
            const SIMD_flt_t delr2x = SIMD_load(tdelx + kcoffset);
            const SIMD_flt_t delr2y = SIMD_load(tdely + kcoffset);
            const SIMD_flt_t delr2z = SIMD_load(tdelz + kcoffset);
            const SIMD_flt_t rsq2 = SIMD_load(trsq + kcoffset);

            const SIMD_flt_t rinvsq2 = SIMD_rcp(rsq2);
            const SIMD_flt_t r2 = SIMD_invsqrt(rinvsq2);
            const SIMD_flt_t rainv2 = SIMD_rcp(r2 - cut);
            const SIMD_flt_t gsrainv2 = sigma_gamma * rainv2;
            const SIMD_flt_t gsrainvsq2 = gsrainv2 * rainv2 / r2;
            const SIMD_flt_t expgsrainv2 = SIMD_exp(gsrainv2);
            const SIMD_flt_t rinv12 = SIMD_rcp(r1 * r2);
            const SIMD_flt_t cs = (delx * delr2x + dely * delr2y +
                              delz * delr2z) * rinv12;
            const SIMD_flt_t delcs = cs - costheta;
            const SIMD_flt_t delcssq = delcs*delcs;

            const SIMD_flt_t facexp = expgsrainv1*expgsrainv2;
            const SIMD_flt_t facrad = lambda_epsilon * facexp * delcssq;
            const SIMD_flt_t frad1 = facrad * gsrainvsq1;
            const SIMD_flt_t frad2 = facrad * gsrainvsq2;
            const SIMD_flt_t facang = lambda_epsilon2 * facexp * delcs;
            const SIMD_flt_t facang12 = rinv12 * facang;
            const SIMD_flt_t csfacang = cs * facang;

            const SIMD_flt_t csfac1 = rinvsq1 * csfacang;
            const SIMD_flt_t fjx = delx * (frad1 + csfac1)-delr2x*facang12;
            const SIMD_flt_t fjy = dely * (frad1 + csfac1)-delr2y*facang12;
            const SIMD_flt_t fjz = delz * (frad1 + csfac1)-delr2z*facang12;

            const SIMD_flt_t csfac2 = rinvsq2 * csfacang;
            SIMD_flt_t fkx = delx * facang12 - delr2x * (frad2 + csfac2);
            SIMD_flt_t fky = dely * facang12 - delr2y * (frad2 + csfac2);
            SIMD_flt_t fkz = delz * facang12 - delr2z * (frad2 + csfac2);

            const SIMD_mask kmask = SIMD_lt(jmask, kk, ejnum);

            SIMD_acc_cache3(kmask, fjx, fjy, fjz, fkx, fky, fkz, fxtmp, fytmp,
                            fztmp, fjxtmp, fjytmp, fjztmp, fxtmp2, fytmp2,
                            fztmp2, fjxtmp2, fjytmp2, fjztmp2,
                            tf + kcoffset * 3, swidth);

            if (EFLAG) {
              SIMD_int k;
              if (eatom) {
                k = SIMD_load(tj + kcoffset);
                k = k << 4;
              }
              SIMD_acc_three(kmask, facrad, eatom, sevdwl, fwtmp, fjtmp,
                             fwtmp2, fjtmp2, k, dforce);
            }
          } // for kk
          if (is_same<flt_t,acc_t>::value == 1)
            SIMD_cache3(tf + coffset * 3, swidth, fjxtmp, fjytmp, fjztmp);
          else
            SIMD_cache3(tf + coffset * 3, swidth, fjxtmp, fjytmp, fjztmp,
                        fjxtmp2, fjytmp2, fjztmp2);

          if (EFLAG) {
            if (eatom) {
              SIMD_int j = SIMD_load(tj + coffset);
              j = j << 4;
              SIMD_jeng_update(jmask, dforce + 3, j, fjtmp);
              if (is_same<flt_t,acc_t>::value == 0)
                SIMD_jeng_update_hi(jmask, dforce + 3, j, fjtmp2);
            }
          }
        } // for jj first loop

        for (int jj = 0; jj < ejnum_max; jj++) {
          const int coffset = jj * swidth;
          const SIMD_mask jmask = jj < ejnum;
          const SIMD_int j = SIMD_load(tj + coffset);
          const SIMD_int joffset = j << 4;

          SIMD_acc_t fjxtmp, fjytmp, fjztmp, fjxtmp2, fjytmp2, fjztmp2;
          int foffset = swidth;
          if (is_same<flt_t,acc_t>::value == 0) foffset = foffset >> 1;
          acc_t *p = tf + coffset * 3;
          fjxtmp = SIMD_load(p);
          if (is_same<flt_t,acc_t>::value == 0) {
            p = p + foffset;
            fjxtmp2 = SIMD_load(p);
          }
          p = p + foffset;
          fjytmp = SIMD_load(p);
          if (is_same<flt_t,acc_t>::value == 0) {
            p = p + foffset;
            fjytmp2 = SIMD_load(p);
          }
          p = p + foffset;
          fjztmp = SIMD_load(p);
          if (is_same<flt_t,acc_t>::value == 0) {
            p = p + foffset;
            fjztmp2 = SIMD_load(p);
          }

          SIMD_conflict_pi_reduce3(jmask, joffset, fjxtmp, fjytmp, fjztmp);
          SIMD_jforce_update(jmask, dforce, joffset, fjxtmp, fjytmp,
                             fjztmp);
          if (is_same<flt_t,acc_t>::value == 0) {
            SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
            SIMD_mask jmask2 = jmask >> 8;
            SIMD_conflict_pi_reduce3(jmask2, joffset2, fjxtmp2, fjytmp2,
                                     fjztmp2);
            SIMD_jforce_update(jmask2, dforce, joffset2, fjxtmp2, fjytmp2,
                               fjztmp2);
          }
        } // for jj second loop

        SIMD_iforce_update(imask, &(f[0].x), goffset, fxtmp, fytmp, fztmp,
                           EFLAG, eatom, fwtmp);
        if (is_same<flt_t,acc_t>::value == 0) {
          imask = imask >> 8;
          SIMD_int goffset2 = _mm512_shuffle_i32x4(goffset, goffset, 238);
          SIMD_iforce_update(imask, &(f[0].x), goffset2, fxtmp2, fytmp2,
                             fztmp2, EFLAG, eatom, fwtmp2);
        }
        if (EFLAG) oevdwl += SIMD_sum(sevdwl);
      } // for ii

      IP_PRE_fdotr_reduce_omp(1, nall, minlocal, nthreads, f_start, f_stride,
                              x, offload, vflag, ov0, ov1, ov2, ov3, ov4, ov5);
    } // end omp

    IP_PRE_fdotr_reduce(1, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      ev_global[0] = oevdwl;
      ev_global[1] = (acc_t)0.0;
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
  } // end offload
  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

#endif

/* ---------------------------------------------------------------------- */

void PairSWIntel::allocate()
{
  PairSW::allocate();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWIntel::init_style()
{
  // there is no support for skipping threebody loops (yet)
  bool tmp_threebody = skip_threebody_flag;
  skip_threebody_flag = false;
  PairSW::init_style();
  skip_threebody_flag = tmp_threebody;

  map[0] = map[1];

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check(true);

  #ifdef _LMP_INTEL_OFFLOAD
  _cop = fix->coprocessor_number();
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    pack_force_const(force_const_single, fix->get_mixed_buffers());
    fix->get_mixed_buffers()->need_tag(1);
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    pack_force_const(force_const_double, fix->get_double_buffers());
    fix->get_double_buffers()->need_tag(1);
  } else {
    pack_force_const(force_const_single, fix->get_single_buffers());
    fix->get_single_buffers()->need_tag(1);
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_noghost())
    error->all(FLERR,"The 'ghost no' option cannot be used with sw/intel.");
  #endif

  #if defined(__INTEL_COMPILER)
  if (__INTEL_COMPILER_BUILD_DATE < 20141023)
    error->all(FLERR, "Intel compiler versions before "
               "15 Update 1 not supported for sw/intel");
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairSWIntel::pack_force_const(ForceConst<flt_t> &fc,
                                   IntelBuffers<flt_t,acc_t> *buffers)
{
  #ifdef LMP_USE_AVXCD
  fix->nbor_pack_width(SIMD_type<flt_t>::width());
  #endif
  fix->three_body_neighbor(1);

  int off_ccache = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop >= 0) off_ccache = 1;
  #endif

  #ifdef LMP_USE_AVXCD
  const int swidth = SIMD_type<flt_t>::width();
  #else
  const int swidth = 1;
  #endif

  buffers->grow_ccache(off_ccache, comm->nthreads, swidth);
  _ccache_stride = buffers->ccache_stride();
  #ifdef LMP_USE_AVXCD
  _ccache_stride3 = buffers->ccache_stride3();
  #endif

  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,memory,_cop);

  // Repeat cutsq calculation because done after call to init_style
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      double cut;
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0))
        cut = init_one(i,j);
      else
        cut = 0.0;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  _onetype = 0;
  _spq = 1;
  int mytypes = 0;
  for (int ii = 0; ii < tp1; ii++) {
    int i = map[ii];
    for (int jj = 0; jj < tp1; jj++) {
      int j = map[jj];
      if (i < 0 || j < 0 || ii == 0 || jj == 0) {
        fc.p2[ii][jj].cutsq = 0;
        fc.p2[ii][jj].cut = 0;
        fc.p2[ii][jj].sigma_gamma = 0;
        fc.p2f[ii][jj].cut = 0;
        fc.p2f[ii][jj].powerp = 0;
        fc.p2f[ii][jj].powerq = 0;
        fc.p2f[ii][jj].sigma = 0;
        fc.p2f2[ii][jj].c1 = 0;
        fc.p2f2[ii][jj].c2 = 0;
        fc.p2f2[ii][jj].c3 = 0;
        fc.p2f2[ii][jj].c4 = 0;
        fc.p2e[ii][jj].c5 = 0;
        fc.p2e[ii][jj].c6 = 0;
      } else {
        int ijparam = elem3param[i][j][j];
        fc.p2[ii][jj].cutsq = params[ijparam].cutsq;
        fc.p2[ii][jj].cut = params[ijparam].cut;
        fc.p2[ii][jj].sigma_gamma = params[ijparam].sigma_gamma;
        fc.p2f[ii][jj].cut = params[ijparam].cut;
        fc.p2f[ii][jj].powerp = -params[ijparam].powerp;
        fc.p2f[ii][jj].powerq = -params[ijparam].powerq;
        fc.p2f[ii][jj].sigma = params[ijparam].sigma;
        fc.p2f2[ii][jj].c1 = params[ijparam].c1;
        fc.p2f2[ii][jj].c2 = params[ijparam].c2;
        fc.p2f2[ii][jj].c3 = params[ijparam].c3;
        fc.p2f2[ii][jj].c4 = params[ijparam].c4;
        fc.p2e[ii][jj].c5 = params[ijparam].c5;
        fc.p2e[ii][jj].c6 = params[ijparam].c6;

        double cutcut = params[ijparam].cut * params[ijparam].cut;
        if (params[ijparam].cutsq >= cutcut)
          fc.p2[ii][jj].cutsq *= 0.98;

        if (params[ijparam].powerp != 4.0 || params[ijparam].powerq != 0.0)
          _spq = 0;
      }

      for (int kk = 0; kk < tp1; kk++) {
        int k = map[kk];
        if (i < 0 || j < 0 || k < 0  || ii == 0 || jj == 0 || kk == 0) {
          fc.p3[ii][jj][kk].costheta = 0;
          fc.p3[ii][jj][kk].lambda_epsilon = 0;
          fc.p3[ii][jj][kk].lambda_epsilon2 = 0;
        } else {
          mytypes++;
          _onetype = ii * tp1 + jj;
          _onetype3 = ii * tp1 * tp1 + jj * tp1 + kk;
          int ijkparam = elem3param[i][j][k];
          fc.p3[ii][jj][kk].costheta = params[ijkparam].costheta;
          fc.p3[ii][jj][kk].lambda_epsilon = params[ijkparam].lambda_epsilon;
          fc.p3[ii][jj][kk].lambda_epsilon2 = params[ijkparam].lambda_epsilon2;
        }
      }
    }
  }
  if (mytypes > 1) _onetype = 0;

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  FC_PACKED0_T *op2 = fc.p2[0];
  FC_PACKED1_T *op2f = fc.p2f[0];
  FC_PACKED1p2_T *op2f2 = fc.p2f2[0];
  FC_PACKED2_T *op2e = fc.p2e[0];
  FC_PACKED3_T *op3 = fc.p3[0][0];
  int tp1sq = tp1 * tp1;
  int tp1cu = tp1sq * tp1;
  #pragma offload_transfer target(mic:_cop) \
    in(op2,op2f,op2f2,op2e: length(tp1sq) alloc_if(0) free_if(0))     \
    in(op3: length(tp1cu) alloc_if(0) free_if(0))
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairSWIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                Memory *memory,
                                                const int cop) {
  if (memory != nullptr) _memory = memory;
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      fc_packed0 *op2 = p2[0];
      fc_packed1 *op2f = p2f[0];
      fc_packed1p2 *op2f2 = p2f2[0];
      fc_packed2 *op2e = p2e[0];
      fc_packed3 *op3 = p3[0][0];

      #ifdef _LMP_INTEL_OFFLOAD
      if (op2 != nullptr && op2f != nullptr && op2f2 != nullptr && op2e != nullptr &&
          op3 != nullptr && _cop >= 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(op2, op2f, op2f2, op2e, op3: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(p2);
      _memory->destroy(p2f);
      _memory->destroy(p2f2);
      _memory->destroy(p2e);
      _memory->destroy(p3);
    }
    if (ntypes > 0) {
      _cop = cop;
      _memory->create(p2,ntypes,ntypes,"fc.p2");
      _memory->create(p2f,ntypes,ntypes,"fc.p2f");
      _memory->create(p2f2,ntypes,ntypes,"fc.p2f2");
      _memory->create(p2e,ntypes,ntypes,"fc.p2e");
      _memory->create(p3,ntypes,ntypes,ntypes,"fc.p3");

      #ifdef _LMP_INTEL_OFFLOAD
      fc_packed0 *op2 = p2[0];
      fc_packed1 *op2f = p2f[0];
      fc_packed1p2 *op2f2 = p2f2[0];
      fc_packed2 *op2e = p2e[0];
      fc_packed3 *op3 = p3[0][0];
      int tp1sq = ntypes * ntypes;
      int tp1cu = tp1sq * ntypes;
      if (op2 != nullptr && op2f != nullptr && op2f2 != nullptr && op2e != nullptr &&
          op3 != nullptr && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(op2,op2f,op2f2,op2e: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(op3: length(tp1cu) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes = ntypes;
}
