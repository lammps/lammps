// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "pair_lj_cut_intel.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>

using namespace LAMMPS_NS;

#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2

/* ---------------------------------------------------------------------- */

PairLJCutIntel::PairLJCutIntel(LAMMPS *lmp) :
  PairLJCut(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJCutIntel::compute(int eflag, int vflag)
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

template <class flt_t, class acc_t>
void PairLJCutIntel::compute(int eflag, int vflag,
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
      buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  int ovflag = 0;
  if (vflag_fdotr) ovflag = 2;
  else if (vflag) ovflag = 1;
  if (_onetype) {
    if (eflag) {
      if (force->newton_pair) {
        eval<1,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (force->newton_pair) {
        eval<1,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  } else {
    if (eflag) {
      if (force->newton_pair) {
        eval<0,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (force->newton_pair) {
        eval<0,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  }
}

template <int ONETYPE, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairLJCutIntel::eval(const int offload, const int vflag,
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

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;
  const flt_t * _noalias const special_lj = fc.special_lj;
  const FC_PACKED1_T * _noalias const ljc12o = fc.ljc12o[0];
  const FC_PACKED2_T * _noalias const lj34 = fc.lj34[0];

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
  int *overflow = fix->get_off_overflow_flag();
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(NEWTON_PAIR, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    if (NEWTON_PAIR == 0 && inum != nlocal)
      memset(f_start, 0, f_stride * sizeof(FORCE_T));

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
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

      flt_t cutsq, lj1, lj2, lj3, lj4, offset;
      if (ONETYPE) {
        cutsq = ljc12o[_onetype].cutsq;
        lj1 = ljc12o[_onetype].lj1;
        lj2 = ljc12o[_onetype].lj2;
        lj3 = lj34[_onetype].lj3;
        lj4 = lj34[_onetype].lj4;
        offset = ljc12o[_onetype].offset;
      }
      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        int itype, ptr_off;
        const FC_PACKED1_T * _noalias ljc12oi;
        const FC_PACKED2_T * _noalias lj34i;
        if (!ONETYPE) {
          itype = x[i].w;
          ptr_off = itype * ntypes;
          ljc12oi = ljc12o + ptr_off;
          lj34i = lj34 + ptr_off;
        }
        const int * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum, offload);

        acc_t fxtmp, fytmp, fztmp, fwtmp;
        acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EFLAG) fwtmp = sevdwl = (acc_t)0;
        if (NEWTON_PAIR == 0)
          if (vflag == VIRIAL_PAIR) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                                   sv0, sv1, sv2, sv3, sv4, sv5)         \
          aligned(jlist,x,ljc12oi,special_lj,f,lj34i:64)
#else
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                               sv0, sv1, sv2, sv3, sv4, sv5)
        #pragma vector aligned
#endif
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          flt_t forcelj, evdwl;
          forcelj = evdwl = (flt_t)0.0;

          int j, jtype, sbindex;
          if (!ONETYPE) {
            sbindex = jlist[jj] >> SBBITS & 3;
            j = jlist[jj] & NEIGHMASK;
          } else
            j = jlist[jj];

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          if (!ONETYPE) {
            jtype = x[j].w;
            cutsq = ljc12oi[jtype].cutsq;
          }
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          #ifdef INTEL_VMASK
          if (rsq < cutsq) {
          #endif
            flt_t factor_lj;
            if (!ONETYPE) factor_lj = special_lj[sbindex];
            flt_t r2inv = 1.0 / rsq;
            flt_t r6inv = r2inv * r2inv * r2inv;
            #ifndef INTEL_VMASK
            if (rsq > cutsq) r6inv = (flt_t)0.0;
            #endif
            if (!ONETYPE) {
              lj1 = ljc12oi[jtype].lj1;
              lj2 = ljc12oi[jtype].lj2;
            }
            forcelj = r6inv * (lj1 * r6inv - lj2);
            flt_t fpair;
            if (!ONETYPE)
              fpair = factor_lj * forcelj * r2inv;
            else
              fpair = forcelj * r2inv;

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
              if (!ONETYPE) {
                lj3 = lj34i[jtype].lj3;
                lj4 = lj34i[jtype].lj4;
                offset = ljc12oi[jtype].offset;
              }
              evdwl = r6inv * (lj3 * r6inv - lj4);
              #ifdef INTEL_VMASK
              evdwl -= offset;
              #else
              if (rsq < cutsq) evdwl -= offset;
              #endif
              if (!ONETYPE) evdwl *= factor_lj;
              sevdwl += evdwl;
              if (eatom) {
                fwtmp += (flt_t)0.5 * evdwl;
                if (NEWTON_PAIR)
                  f[j].w += (flt_t)0.5 * evdwl;
              }
            }

            if (NEWTON_PAIR == 0)
              IP_PRE_ev_tally_nborv(vflag, delx, dely, delz, fpx, fpy, fpz);
          #ifdef INTEL_VMASK
          } // if rsq
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

        IP_PRE_ev_tally_atom(NEWTON_PAIR, EFLAG, vflag, f, fwtmp);
      } // for ii

      IP_PRE_fdotr_reduce_omp(NEWTON_PAIR, nall, minlocal, nthreads, f_start,
                              f_stride, x, offload, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
    } // end omp

    IP_PRE_fdotr_reduce(NEWTON_PAIR, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      if (NEWTON_PAIR == 0) {
        oevdwl *= (acc_t)0.5;
        ov0 *= (acc_t)0.5;
        ov1 *= (acc_t)0.5;
        ov2 *= (acc_t)0.5;
        ov3 *= (acc_t)0.5;
        ov4 *= (acc_t)0.5;
        ov5 *= (acc_t)0.5;
      }
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

/* ---------------------------------------------------------------------- */

void PairLJCutIntel::init_style()
{
  PairLJCut::init_style();
  auto request = neighbor->find_request(this);

  if (force->newton_pair == 0) {
    request->half = 0;
    request->full = 1;
  }
  request->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_balance() != 0.0)
    error->all(FLERR,
          "Offload for lj/cut/intel is not yet available. Set balance to 0.");
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairLJCutIntel::pack_force_const(ForceConst<flt_t> &fc,
                                      IntelBuffers<flt_t,acc_t> *buffers)
{
  _onetype = 0;

  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,memory,_cop);

  // Repeat cutsq calculation because done after call to init_style
  int mytypes = 0;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      double cut;
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        mytypes++;
        _onetype = i * tp1 + j;
      } else {
        cut = 0.0;
      }
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }
  if (mytypes > 1 || atom->molecular) _onetype = 0;

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.ljc12o[i][j].lj1 = lj1[i][j];
      fc.ljc12o[i][j].lj2 = lj2[i][j];
      fc.lj34[i][j].lj3 = lj3[i][j];
      fc.lj34[i][j].lj4 = lj4[i][j];
      fc.ljc12o[i][j].cutsq = cutsq[i][j];
      fc.ljc12o[i][j].offset = offset[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairLJCutIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                   Memory *memory,
                                                   const int cop) {
  if (memory != nullptr) _memory = memory;
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      _memory->destroy(ljc12o);
      _memory->destroy(lj34);
    }
    if (ntypes > 0) {
      _cop = cop;
      _memory->create(ljc12o,ntypes,ntypes,"fc.c12o");
      _memory->create(lj34,ntypes,ntypes,"fc.lj34");
    }
  }
  _ntypes = ntypes;
}
