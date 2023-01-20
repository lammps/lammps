// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "pair_lj_charmm_coul_charmm_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define LJ_T typename IntelBuffers<flt_t,flt_t>::vec4_t

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulCharmmIntel::PairLJCharmmCoulCharmmIntel(LAMMPS *lmp) :
  PairLJCharmmCoulCharmm(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulCharmmIntel::compute(int eflag, int vflag)
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
void PairLJCharmmCoulCharmmIntel::compute(int eflag, int vflag,
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
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal+atom->nghost,
                                packthreads, sizeof(ATOM_T));
      buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  // -------------------- Regular version
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
void PairLJCharmmCoulCharmmIntel::eval(const int offload, const int vflag,
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
  const flt_t inv_denom_lj = (flt_t)1.0/denom_lj;
  const flt_t inv_denom_coul = (flt_t)1.0/denom_coul;

  const flt_t * _noalias const cutsq = fc.cutsq[0];
  const LJ_T * _noalias const lj = fc.lj[0];
  const flt_t cut_ljsq = fc.cut_ljsq;
  const flt_t cut_lj_innersq = fc.cut_lj_innersq;
  const flt_t cut_coul_innersq = fc.cut_coul_innersq;
  const flt_t cut_coulsq = fc.cut_coulsq;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  flt_t * _noalias const ccachex = buffers->get_ccachex();
  flt_t * _noalias const ccachey = buffers->get_ccachey();
  flt_t * _noalias const ccachez = buffers->get_ccachez();
  flt_t * _noalias const ccachew = buffers->get_ccachew();
  int * _noalias const ccachei = buffers->get_ccachei();
  int * _noalias const ccachej = buffers->get_ccachej();
  const int ccache_stride = _ccache_stride;

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

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if (offload) \
    in(special_lj,special_coul:length(0) alloc_if(0) free_if(0)) \
    in(cutsq,lj:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(q:length(q_size) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(ccachex,ccachey,ccachez,ccachew:length(0) alloc_if(0) free_if(0)) \
    in(ccachei,ccachej:length(0) alloc_if(0) free_if(0)) \
    in(ccache_stride,nthreads,qqrd2e,inum,nall,ntypes,cut_coulsq) \
    in(vflag,eatom,f_stride,separate_flag,offload) \
    in(astart,cut_ljsq,cut_lj_innersq,nlocal,inv_denom_lj,minlocal) \
    in(inv_denom_coul,cut_coul_innersq) \
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
      flt_t cutboth = cut_coulsq;

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;

      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        const int itype = x[i].w;

        const int ptr_off = itype * ntypes;
        const flt_t * _noalias const cutsqi = cutsq + ptr_off;
        const LJ_T * _noalias const lji = lj + ptr_off;

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
          if (vflag == VIRIAL_PAIR) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

        int ej = 0;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          const int j = jlist[jj] & NEIGHMASK;
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cut_coulsq) {
            trsq[ej]=rsq;
            tdelx[ej]=delx;
            tdely[ej]=dely;
            tdelz[ej]=delz;
            tjtype[ej]=x[j].w;
            tj[ej]=jlist[jj];
            ej++;
          }
        }

        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                                   secoul, sv0, sv1, sv2, sv3, sv4, sv5)
#else
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                               secoul, sv0, sv1, sv2, sv3, sv4, sv5)
#endif
        #pragma vector aligned
        #endif
        for (int jj = 0; jj < ej; jj++) {
          flt_t forcecoul, forcelj, evdwl;
          forcecoul = forcelj = evdwl = (flt_t)0.0;

          const int j = tj[jj] & NEIGHMASK;
          const int sbindex = tj[jj] >> SBBITS & 3;
          const flt_t rsq = trsq[jj];
          const flt_t r2inv = (flt_t)1.0 / rsq;
          const flt_t r_inv = (flt_t)1.0 / sqrt(rsq);
          forcecoul = qqrd2e * qtmp * q[j] * r_inv;
          if (rsq > cut_coul_innersq) {
            const flt_t ccr = cut_coulsq - rsq;
            const flt_t switch1 = ccr * ccr * inv_denom_coul *
              (cut_coulsq + (flt_t)2.0 * rsq - (flt_t)3.0 * cut_coul_innersq);
            forcecoul *= switch1;
          }

          #ifdef INTEL_VMASK
          if (rsq < cut_ljsq) {
          #endif
            const int jtype = IP_PRE_dword_index(tjtype[jj]);
            flt_t r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lji[jtype].x * r6inv - lji[jtype].y);
            if (EFLAG) evdwl = r6inv*(lji[jtype].z * r6inv - lji[jtype].w);

            #ifdef INTEL_VMASK
            if (rsq > cut_lj_innersq) {
            #endif
              const flt_t drsq = cut_ljsq - rsq;
              const flt_t cut2 = (rsq - cut_lj_innersq) * drsq;
              const flt_t switch1 = drsq * (drsq * drsq + (flt_t)3.0 * cut2) *
                  inv_denom_lj;
              const flt_t switch2 = (flt_t)12.0 * rsq * cut2 * inv_denom_lj;
              if (EFLAG) {
                #ifndef INTEL_VMASK
                if (rsq > cut_lj_innersq) {
                #endif
                  forcelj = forcelj * switch1 + evdwl * switch2;
                  evdwl *= switch1;
                #ifndef INTEL_VMASK
                }
                #endif
              } else {
                const flt_t philj = r6inv * (lji[jtype].z*r6inv -
                    lji[jtype].w);
                #ifndef INTEL_VMASK
                if (rsq > cut_lj_innersq)
                #endif
                  forcelj =  forcelj * switch1 + philj * switch2;
              }
            #ifdef INTEL_VMASK
            }
            #endif

          #ifdef INTEL_VMASK
          }
          #else
          if (rsq > cut_ljsq) { forcelj = (flt_t)0.0; evdwl = (flt_t)0.0; }
          #endif
          if (sbindex) {
            const flt_t factor_coul = special_coul[sbindex];
            forcecoul *= factor_coul;
            const flt_t factor_lj = special_lj[sbindex];
            forcelj *= factor_lj;
            if (EFLAG) evdwl *= factor_lj;
          }

          const flt_t fpair = (forcecoul + forcelj) * r2inv;
          const flt_t fpx = fpair * tdelx[jj];
          fxtmp += fpx;
          if (NEWTON_PAIR) f[j].x -= fpx;
          const flt_t fpy = fpair * tdely[jj];
          fytmp += fpy;
          if (NEWTON_PAIR) f[j].y -= fpy;
          const flt_t fpz = fpair * tdelz[jj];
          fztmp += fpz;
          if (NEWTON_PAIR) f[j].z -= fpz;

          if (EFLAG) {
            sevdwl += evdwl;
            secoul += forcecoul;
            if (eatom) {
              fwtmp += (flt_t)0.5 * evdwl + (flt_t)0.5 * forcecoul;
              if (NEWTON_PAIR)
                f[j].w += (flt_t)0.5 * evdwl + (flt_t)0.5 * forcecoul;
            }
          }
          if (NEWTON_PAIR == 0)
            IP_PRE_ev_tally_nborv(vflag, tdelx[jj], tdely[jj], tdelz[jj],
                                  fpx, fpy, fpz);
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

void PairLJCharmmCoulCharmmIntel::init_style()
{
  PairLJCharmmCoulCharmm::init_style();
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
void PairLJCharmmCoulCharmmIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int off_ccache = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop >= 0) off_ccache = 1;
  #endif
  buffers->grow_ccache(off_ccache, comm->nthreads, 1);
  _ccache_stride = buffers->ccache_stride();

  int tp1 = atom->ntypes + 1;

  fc.set_ntypes(tp1, memory, _cop);

  // Repeat cutsq calculation because done after call to init_style
  if (cut_lj > cut_coul)
    error->all(FLERR,
         "Intel variant of lj/charmm/coul/long expects lj cutoff<=coulombic");
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

  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  fc.cut_coulsq = cut_coulsq;
  fc.cut_ljsq = cut_ljsq;
  fc.cut_coul_innersq = cut_coul_innersq;
  fc.cut_lj_innersq = cut_lj_innersq;

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_coul[i] = force->special_coul[i];
    fc.special_coul[0] = 1.0;
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.lj[i][j].x = lj1[i][j];
      fc.lj[i][j].y = lj2[i][j];
      fc.lj[i][j].z = lj3[i][j];
      fc.lj[i][j].w = lj4[i][j];
      fc.cutsq[i][j] = cutsq[i][j];
    }
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  flt_t * special_lj = fc.special_lj;
  flt_t * special_coul = fc.special_coul;
  flt_t * cutsq = fc.cutsq[0];
  LJ_T * lj = fc.lj[0];
  int tp1sq = tp1 * tp1;
  #pragma offload_transfer target(mic:_cop) \
    in(special_lj, special_coul: length(4) alloc_if(0) free_if(0)) \
    in(cutsq,lj: length(tp1sq) alloc_if(0) free_if(0))
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairLJCharmmCoulCharmmIntel::ForceConst<flt_t>::set_ntypes(
  const int ntypes, Memory *memory, const int cop) {
  if (memory != nullptr) _memory = memory;
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      flt_t * ocutsq = cutsq[0];
      typename IntelBuffers<flt_t,flt_t>::vec4_t * olj = lj[0];
      if (ospecial_lj != nullptr && ocutsq != nullptr && olj != nullptr &&
          ospecial_coul != nullptr && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj, ospecial_coul: alloc_if(0) free_if(1)) \
          nocopy(ocutsq, olj: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(cutsq);
      _memory->destroy(lj);
    }
    if (ntypes > 0) {
      _cop = cop;
      _memory->create(cutsq,ntypes,ntypes,"fc.cutsq");
      _memory->create(lj,ntypes,ntypes,"fc.lj");

      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      flt_t * ocutsq = cutsq[0];
      typename IntelBuffers<flt_t,flt_t>::vec4_t * olj = lj[0];
      int tp1sq = ntypes*ntypes;
      if (ospecial_lj != nullptr && ocutsq != nullptr && olj != nullptr &&
          ospecial_coul != nullptr && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj: length(4) alloc_if(1) free_if(0)) \
          nocopy(ospecial_coul: length(4) alloc_if(1) free_if(0)) \
          nocopy(ocutsq,olj: length(tp1sq) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes=ntypes;
}
