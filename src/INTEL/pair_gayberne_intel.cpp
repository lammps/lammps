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

#include "pair_gayberne_intel.h"
#include "math_extra_intel.h"

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif
#include <cmath>
#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cstring>

using namespace LAMMPS_NS;

#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2
#define FC_PACKED3_T typename ForceConst<flt_t>::fc_packed3

/* ---------------------------------------------------------------------- */

PairGayBerneIntel::PairGayBerneIntel(LAMMPS *lmp) :
  PairGayBerne(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairGayBerneIntel::compute(int eflag, int vflag)
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
void PairGayBerneIntel::compute(int eflag, int vflag,
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
  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int host_start = fix->host_start_pair();
  const int offload_end = fix->offload_end_pair();
  const int ago = neighbor->ago;

  if (fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);
    const AtomVecEllipsoid::Bonus * const bonus = avec->bonus;
    const int * const ellipsoid = atom->ellipsoid;
    QUAT_T * _noalias const quat = buffers->get_quat();

    int packthreads;
    if (nthreads > INTEL_HTHREADS) packthreads = nthreads;
    else packthreads = 1;
    #if defined(_OPENMP)
    #pragma omp parallel if (packthreads > 1)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, nall, packthreads,
                                sizeof(ATOM_T));
      if (ago != 0) buffers->thr_pack(ifrom,ito,ago);

      for (int i = ifrom; i < ito; i++) {
        int qi = ellipsoid[i];
        if (qi > -1) {
          quat[i].w = bonus[qi].quat[0];
          quat[i].i = bonus[qi].quat[1];
          quat[i].j = bonus[qi].quat[2];
          quat[i].k = bonus[qi].quat[3];
        }
      }
    }
    quat[nall].w = (flt_t)1.0;
    quat[nall].i = (flt_t)0.0;
    quat[nall].j = (flt_t)0.0;
    quat[nall].k = (flt_t)0.0;
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

template <int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairGayBerneIntel::eval(const int offload, const int vflag,
                             IntelBuffers<flt_t,acc_t> *buffers,
                             const ForceConst<flt_t> &fc,
                             const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  ATOM_T * _noalias const x = buffers->get_x(offload);
  QUAT_T * _noalias const quat = buffers->get_quat(offload);
  const AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  const int *ellipsoid = atom->ellipsoid;

  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->separate_buffers()) {
    fix->start_watch(TIME_PACK);
    if (offload) {
      #pragma omp parallel
      {
        int ifrom, ito, tid;
        int nthreads = comm->nthreads;
        IP_PRE_omp_range_id_align(ifrom, ito, tid, nlocal,
                                  nthreads, sizeof(ATOM_T));
        if (ago != 0) buffers->thr_pack_cop(ifrom, ito, 0);
        for (int i = ifrom; i < ito; i++) {
          int qi = ellipsoid[i];
          if (qi > -1) {
            quat[i].w = bonus[qi].quat[0];
            quat[i].i = bonus[qi].quat[1];
            quat[i].j = bonus[qi].quat[2];
            quat[i].k = bonus[qi].quat[3];
          }
        }
        int nghost = nall - nlocal;
        if (nghost) {
          IP_PRE_omp_range_align(ifrom, ito, tid, nall - nlocal,
                                 nthreads, sizeof(ATOM_T));
          int offset = 0;
          ifrom += nlocal;
          ito += nlocal;
          if (ago != 0) {
            offset = fix->offload_min_ghost() - nlocal;
            buffers->thr_pack_cop(ifrom, ito, offset, ago == 1);
          }
          for (int i = ifrom; i < ito; i++) {
            int qi = ellipsoid[i + offset];
            if (qi > -1) {
              quat[i].w = bonus[qi].quat[0];
              quat[i].i = bonus[qi].quat[1];
              quat[i].j = bonus[qi].quat[2];
              quat[i].k = bonus[qi].quat[3];
            }
          }
        }
      }
    } else {
      if (ago != 0) buffers->thr_pack_host(fix->host_min_local(), nlocal, 0);
      for (int i = fix->host_min_local(); i < nlocal; i++) {
        int qi = ellipsoid[i];
        if (qi > -1) {
          quat[i].w = bonus[qi].quat[0];
          quat[i].i = bonus[qi].quat[1];
          quat[i].j = bonus[qi].quat[2];
          quat[i].k = bonus[qi].quat[3];
        }
      }
      int offset = fix->host_min_ghost() - nlocal;
      if (ago != 0) buffers->thr_pack_host(nlocal, nall, offset);
      for (int i = nlocal; i < nall; i++) {
        int qi = ellipsoid[i + offset];
        if (qi > -1) {
          quat[i].w = bonus[qi].quat[0];
          quat[i].i = bonus[qi].quat[1];
          quat[i].j = bonus[qi].quat[2];
          quat[i].k = bonus[qi].quat[3];
        }
      }
    }
    fix->stop_watch(TIME_PACK);
  }
  #endif

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;  // NOLINT
  const flt_t * _noalias const special_lj = fc.special_lj;

  const FC_PACKED1_T * _noalias const ijc = fc.ijc[0];
  const FC_PACKED2_T * _noalias const lj34 = fc.lj34[0];
  const FC_PACKED3_T * _noalias const ic = fc.ic;
  const flt_t mu = fc.mu;
  const flt_t gamma = fc.gamma;
  const flt_t upsilon = fc.upsilon;

  flt_t * const rsq_formi = fc.rsq_form[0];
  flt_t * const delx_formi = fc.delx_form[0];
  flt_t * const dely_formi = fc.dely_form[0];
  flt_t * const delz_formi = fc.delz_form[0];
  int * const jtype_formi = fc.jtype_form[0];
  int * const jlist_formi = fc.jlist_form[0];

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
  const int max_nbors = _max_nbors;
  const int nthreads = tc;

  int pad = 1;
  if (offload) {
    if (INTEL_MIC_NBOR_PAD > 1)
      pad = INTEL_MIC_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  } else {
    if (INTEL_NBOR_PAD > 1)
      pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  }
  const int pad_width = pad;

  #ifdef _LMP_INTEL_OFFLOAD
  int *overflow = fix->get_off_overflow_flag();
  double *timer_compute = fix->off_watch_pair();

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if (offload) \
    in(special_lj:length(0) alloc_if(0) free_if(0)) \
    in(ijc,lj34,ic:length(0) alloc_if(0) free_if(0)) \
    in(rsq_formi, delx_formi, dely_formi: length(0) alloc_if(0) free_if(0)) \
    in(delz_formi, jtype_formi, jlist_formi: length(0) alloc_if(0) free_if(0))\
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(quat:length(nall+1) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(nthreads,inum,nall,ntypes,vflag,eatom,minlocal,separate_flag) \
    in(astart,nlocal,f_stride,max_nbors,mu,gamma,upsilon,offload,pad_width) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute=MIC_Wtime();
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    if (separate_flag) {
      if (separate_flag < 3) {
        int all_local = nlocal;
        int ghost_min = overflow[LMP_GHOST_MIN];
        nlocal = overflow[LMP_LOCAL_MAX] + 1;
        int nghost = overflow[LMP_GHOST_MAX] + 1 - ghost_min;
        if (nghost < 0) nghost = 0;
        nall = nlocal + nghost;
        separate_flag--;
        int flength;
        if (NEWTON_PAIR) flength = nall;
        else flength = nlocal;
        IP_PRE_get_stride(f_stride, flength, sizeof(FORCE_T),
                             separate_flag);
        if (nghost) {
          if (nlocal < all_local || ghost_min > all_local) {
            memmove(x + nlocal, x + ghost_min,
                    (nall - nlocal) * sizeof(ATOM_T));
            memmove(quat + nlocal, quat + ghost_min,
                    (nall - nlocal) * sizeof(QUAT_T));
          }
        }
      }
      x[nall].x = (flt_t)INTEL_BIGP;
      x[nall].y = (flt_t)INTEL_BIGP;
      x[nall].z = (flt_t)INTEL_BIGP;
      x[nall].w = 1;
      quat[nall].w = (flt_t)1.0;
      quat[nall].i = (flt_t)0.0;
      quat[nall].j = (flt_t)0.0;
      quat[nall].k = (flt_t)0.0;
    }
    #endif

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
    if (NEWTON_PAIR == 0) f_start[1].w = 0;
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
      if (NEWTON_PAIR) foff = tid * f_stride - minlocal * 2;
      else foff = minlocal*-2;
      FORCE_T * _noalias const f = f_start + foff;
      if (NEWTON_PAIR) memset(f + minlocal * 2, 0, f_stride * sizeof(FORCE_T));

      flt_t * _noalias const rsq_form = rsq_formi + tid * max_nbors;
      flt_t * _noalias const delx_form = delx_formi + tid * max_nbors;
      flt_t * _noalias const dely_form = dely_formi + tid * max_nbors;
      flt_t * _noalias const delz_form = delz_formi + tid * max_nbors;
      int * _noalias const jtype_form = jtype_formi + tid * max_nbors;
      int * _noalias const jlist_form = jlist_formi + tid * max_nbors;

      int ierror = 0;
      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        const int itype = x[i].w;
        const int ptr_off = itype * ntypes;
        const FC_PACKED1_T * _noalias const ijci = ijc + ptr_off;
        const FC_PACKED2_T * _noalias const lj34i = lj34 + ptr_off;

        const int * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum, offload);

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;

        flt_t a1_0, a1_1, a1_2, a1_3, a1_4, a1_5, a1_6, a1_7, a1_8;
        flt_t b1_0, b1_1, b1_2, b1_3, b1_4, b1_5, b1_6, b1_7, b1_8;
        flt_t g1_0, g1_1, g1_2, g1_3, g1_4, g1_5, g1_6, g1_7, g1_8;

        if (ijci[itype].form == ELLIPSE_ELLIPSE) {
          flt_t temp_0,temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7,temp_8;
          ME_quat_to_mat_trans(quat[i],a1);
          ME_diag_times3(ic[itype].well,a1,temp);
          ME_transpose_times3(a1,temp,b1);
          ME_diag_times3(ic[itype].shape2,a1,temp);
          ME_transpose_times3(a1,temp,g1);
        }

        acc_t fxtmp, fytmp, fztmp, fwtmp, t1tmp, t2tmp, t3tmp;
        acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;
        fxtmp = fytmp = fztmp = t1tmp = t2tmp = t3tmp = (acc_t)0.0;

        if (EFLAG) fwtmp = sevdwl = (acc_t)0.0;
        if (NEWTON_PAIR == 0)
          if (vflag == VIRIAL_PAIR) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;

        bool multiple_forms = false;
        int packed_j = 0;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          int jm = jlist[jj];
          int j = jm & NEIGHMASK;
          const int jtype = IP_PRE_dword_index(x[j].w);

          if (ijci[jtype].form == ELLIPSE_ELLIPSE) {
            flt_t delx = x[j].x-xtmp;
            flt_t dely = x[j].y-ytmp;
            flt_t delz = x[j].z-ztmp;
            flt_t rsq = delx * delx + dely * dely + delz * delz;

            if (rsq < ijci[jtype].cutsq) {
              rsq_form[packed_j] = rsq;
              delx_form[packed_j] = delx;
              dely_form[packed_j] = dely;
              delz_form[packed_j] = delz;
              jtype_form[packed_j] = jtype;
              jlist_form[packed_j] = jm;
              packed_j++;
            }
          } else
            multiple_forms = true;
        }
        int packed_end = packed_j;
        IP_PRE_neighbor_pad(packed_end, offload);
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min=1, max=15, avg=8
        #endif
        for ( ; packed_j < packed_end; packed_j++)
          jlist_form[packed_j] = nall;

        // -------------------------------------------------------------

        #ifdef INTEL_V512
        __assume(packed_j % INTEL_VECTOR_WIDTH == 0);
        __assume(packed_j % 8 == 0);
        __assume(packed_j % INTEL_MIC_VECTOR_WIDTH == 0);
        #endif
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd reduction(+:fxtmp,fytmp,fztmp,fwtmp,t1tmp,t2tmp, \
                                   t3tmp,sevdwl,sv0,sv1,sv2,sv3,sv4,sv5)
#else
        #pragma simd reduction(+:fxtmp,fytmp,fztmp,fwtmp,t1tmp,t2tmp, \
                               t3tmp,sevdwl,sv0,sv1,sv2,sv3,sv4,sv5)
#endif
        #pragma vector aligned
        #endif
        for (int jj = 0; jj < packed_j; jj++) {
          flt_t a2_0, a2_1, a2_2, a2_3, a2_4, a2_5, a2_6, a2_7, a2_8;
          flt_t b2_0, b2_1, b2_2, b2_3, b2_4, b2_5, b2_6, b2_7, b2_8;
          flt_t g2_0, g2_1, g2_2, g2_3, g2_4, g2_5, g2_6, g2_7, g2_8;
          flt_t temp_0,temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7,temp_8;
          flt_t fforce_0, fforce_1, fforce_2, ttor_0, ttor_1, ttor_2;
          flt_t rtor_0, rtor_1, rtor_2;

          const int sbindex = jlist_form[jj] >> SBBITS & 3;
          const int j = jlist_form[jj] & NEIGHMASK;
          flt_t factor_lj = special_lj[sbindex];
          const int jtype = IP_PRE_dword_index(jtype_form[jj]);
          const flt_t sigma = ijci[jtype].sigma;
          const flt_t epsilon = ijci[jtype].epsilon;
          const flt_t shape2_0 = ic[jtype].shape2[0];
          const flt_t shape2_1 = ic[jtype].shape2[1];
          const flt_t shape2_2 = ic[jtype].shape2[2];
          flt_t one_eng, evdwl;

          ME_quat_to_mat_trans(quat[j], a2);
          ME_diag_times3(ic[jtype].well, a2, temp);
          ME_transpose_times3(a2, temp, b2);
          ME_diag_times3a(shape2, a2, temp);
          ME_transpose_times3(a2, temp, g2);

          flt_t tempv_0, tempv_1, tempv_2, tempv2_0, tempv2_1, tempv2_2;
          flt_t temp1, temp2, temp3;

          flt_t r12hat_0, r12hat_1, r12hat_2;
          ME_normalize3(delx_form[jj], dely_form[jj], delz_form[jj], r12hat);
          flt_t r = std::sqrt(rsq_form[jj]);

          // compute distance of closest approach

          flt_t g12_0, g12_1, g12_2, g12_3, g12_4, g12_5, g12_6, g12_7, g12_8;
          ME_plus3(g1, g2, g12);
          flt_t kappa_0, kappa_1, kappa_2;
          ME_mldivide3(g12, delx_form[jj], dely_form[jj], delz_form[jj],
                       kappa, ierror);

          // tempv = G12^-1*r12hat

          flt_t inv_r = (flt_t)1.0 / r;
          tempv_0 = kappa_0 * inv_r;
          tempv_1 = kappa_1 * inv_r;
          tempv_2 = kappa_2 * inv_r;
          flt_t sigma12 = ME_dot3(r12hat, tempv);
          sigma12 = std::pow((flt_t)0.5 * sigma12,(flt_t) - 0.5);
          flt_t h12 = r - sigma12;

          // energy
          // compute u_r

          flt_t varrho = sigma / (h12 + gamma * sigma);
          flt_t varrho6 = std::pow(varrho, (flt_t)6.0);
          flt_t varrho12 = varrho6 * varrho6;
          flt_t u_r = (flt_t)4.0 * epsilon * (varrho12 - varrho6);

          // compute eta_12

          flt_t eta = (flt_t)2.0 * ijci[jtype].lshape;
          flt_t det_g12 = ME_det3(g12);
          eta = std::pow(eta / det_g12, upsilon);

          // compute chi_12

          flt_t b12_0, b12_1, b12_2, b12_3, b12_4, b12_5, b12_6, b12_7, b12_8;
          flt_t iota_0, iota_1, iota_2;
          ME_plus3(b1, b2, b12);
          ME_mldivide3(b12, delx_form[jj], dely_form[jj], delz_form[jj],
                       iota, ierror);

          // tempv = G12^-1*r12hat

          tempv_0 = iota_0 * inv_r;
          tempv_1 = iota_1 * inv_r;
          tempv_2 = iota_2 * inv_r;
          flt_t chi = ME_dot3(r12hat, tempv);
          chi = std::pow(chi * (flt_t)2.0, mu);

          // force
          // compute dUr/dr

          temp1 = ((flt_t)2.0 * varrho12 * varrho - varrho6 * varrho) /
            sigma;
          temp1 = temp1 * (flt_t)24.0 * epsilon;
          flt_t u_slj = temp1 * std::pow(sigma12, (flt_t)3.0) * (flt_t)0.5;
          flt_t dUr_0, dUr_1, dUr_2;
          temp2 = ME_dot3(kappa, r12hat);
          flt_t uslj_rsq = u_slj / rsq_form[jj];
          dUr_0 = temp1 * r12hat_0 + uslj_rsq * (kappa_0 - temp2 * r12hat_0);
          dUr_1 = temp1 * r12hat_1 + uslj_rsq * (kappa_1 - temp2 * r12hat_1);
          dUr_2 = temp1 * r12hat_2 + uslj_rsq * (kappa_2 - temp2 * r12hat_2);

          // compute dChi_12/dr

          flt_t dchi_0, dchi_1, dchi_2;
          temp1 = ME_dot3(iota, r12hat);
          temp2 = (flt_t)-4.0 / rsq_form[jj] * mu *
            std::pow(chi, (mu - (flt_t)1.0) / mu);
          dchi_0 = temp2 * (iota_0 - temp1 * r12hat_0);
          dchi_1 = temp2 * (iota_1 - temp1 * r12hat_1);
          dchi_2 = temp2 * (iota_2 - temp1 * r12hat_2);

          temp1 = -eta * u_r;
          temp3 = eta * chi;
          fforce_0 = temp1 * dchi_0 - temp3 * dUr_0;
          fforce_1 = temp1 * dchi_1 - temp3 * dUr_1;
          fforce_2 = temp1 * dchi_2 - temp3 * dUr_2;

          // torque for particle 1 and 2
          // compute dUr

          tempv_0 = -uslj_rsq * kappa_0;
          tempv_1 = -uslj_rsq * kappa_1;
          tempv_2 = -uslj_rsq * kappa_2;
          ME_vecmat(kappa, g1, tempv2);
          ME_cross3(tempv, tempv2, dUr);
          flt_t dUr2_0, dUr2_1, dUr2_2;

          if (NEWTON_PAIR) {
            ME_vecmat(kappa, g2, tempv2);
            ME_cross3(tempv, tempv2, dUr2);
          }

          // compute d_chi

          ME_vecmat(iota, b1, tempv);
          ME_cross3(tempv, iota, dchi);
          dchi_0 *= temp2;
          dchi_1 *= temp2;
          dchi_2 *= temp2;
          flt_t dchi2_0, dchi2_1, dchi2_2;

          if (NEWTON_PAIR) {
            ME_vecmat(iota, b2, tempv);
            ME_cross3(tempv, iota, dchi2);
            dchi2_0 *= temp2;
            dchi2_1 *= temp2;
            dchi2_2 *= temp2;
          }

          // compute d_eta

          flt_t deta_0, deta_1, deta_2;
          deta_0 = deta_1 = deta_2 = (flt_t)0.0;
          ME_compute_eta_torque(g12, a1, shape2, temp);
          temp1 = -eta * upsilon;

          tempv_0 = temp1 * temp_0;
          tempv_1 = temp1 * temp_1;
          tempv_2 = temp1 * temp_2;
          ME_mv0_cross3(a1, tempv, tempv2);
          deta_0 += tempv2_0;
          deta_1 += tempv2_1;
          deta_2 += tempv2_2;

          tempv_0 = temp1 * temp_3;
          tempv_1 = temp1 * temp_4;
          tempv_2 = temp1 * temp_5;
          ME_mv1_cross3(a1, tempv, tempv2);
          deta_0 += tempv2_0;
          deta_1 += tempv2_1;
          deta_2 += tempv2_2;

          tempv_0 = temp1 * temp_6;
          tempv_1 = temp1 * temp_7;
          tempv_2 = temp1 * temp_8;
          ME_mv2_cross3(a1, tempv, tempv2);
          deta_0 += tempv2_0;
          deta_1 += tempv2_1;
          deta_2 += tempv2_2;

          // compute d_eta for particle 2

          flt_t deta2_0, deta2_1, deta2_2;
          if (NEWTON_PAIR) {
            deta2_0 = deta2_1 = deta2_2 = (flt_t)0.0;
            ME_compute_eta_torque(g12, a2, shape2, temp);

            tempv_0 = temp1 * temp_0;
            tempv_1 = temp1 * temp_1;
            tempv_2 = temp1 * temp_2;
            ME_mv0_cross3(a2, tempv, tempv2);
            deta2_0 += tempv2_0;
            deta2_1 += tempv2_1;
            deta2_2 += tempv2_2;

            tempv_0 = temp1 * temp_3;
            tempv_1 = temp1 * temp_4;
            tempv_2 = temp1 * temp_5;
            ME_mv1_cross3(a2, tempv, tempv2);
            deta2_0 += tempv2_0;
            deta2_1 += tempv2_1;
            deta2_2 += tempv2_2;

            tempv_0 = temp1 * temp_6;
            tempv_1 = temp1 * temp_7;
            tempv_2 = temp1 * temp_8;
            ME_mv2_cross3(a2, tempv, tempv2);
            deta2_0 += tempv2_0;
            deta2_1 += tempv2_1;
            deta2_2 += tempv2_2;
          }

          // torque

          temp1 = u_r * eta;
          temp2 = u_r * chi;
          temp3 = chi * eta;

          ttor_0 = (temp1 * dchi_0 + temp2 * deta_0 + temp3 * dUr_0) *
            (flt_t)-1.0;
          ttor_1 = (temp1 * dchi_1 + temp2 * deta_1 + temp3 * dUr_1) *
            (flt_t)-1.0;
          ttor_2 = (temp1 * dchi_2 + temp2 * deta_2 + temp3 * dUr_2) *
            (flt_t)-1.0;

          if (NEWTON_PAIR) {
            rtor_0 = (temp1 * dchi2_0 + temp2 * deta2_0 + temp3 * dUr2_0) *
              (flt_t)-1.0;
            rtor_1 = (temp1 * dchi2_1 + temp2 * deta2_1 + temp3 * dUr2_1) *
              (flt_t)-1.0;
            rtor_2 = (temp1 * dchi2_2 + temp2 * deta2_2 + temp3 * dUr2_2) *
              (flt_t)-1.0;
          }

          one_eng = temp1 * chi;
          #ifndef INTEL_VMASK
          if (jlist_form[jj] == nall) {
            one_eng = (flt_t)0.0;
            fforce_0 = 0.0;
            fforce_1 = 0.0;
            fforce_2 = 0.0;
            ttor_0 = 0.0;
            ttor_1 = 0.0;
            ttor_2 = 0.0;
            rtor_0 = 0.0;
            rtor_1 = 0.0;
            rtor_2 = 0.0;
          }
          #endif

          fforce_0 *= factor_lj;
          fforce_1 *= factor_lj;
          fforce_2 *= factor_lj;
          ttor_0 *= factor_lj;
          ttor_1 *= factor_lj;
          ttor_2 *= factor_lj;

          #ifdef INTEL_VMASK
          if (jlist_form[jj] < nall) {
          #endif
            fxtmp += fforce_0;
            fytmp += fforce_1;
            fztmp += fforce_2;
            t1tmp += ttor_0;
            t2tmp += ttor_1;
            t3tmp += ttor_2;

            if (NEWTON_PAIR) {
              rtor_0 *= factor_lj;
              rtor_1 *= factor_lj;
              rtor_2 *= factor_lj;
              int jp = j * 2;
              f[jp].x -= fforce_0;
              f[jp].y -= fforce_1;
              f[jp].z -= fforce_2;
              jp++;
              f[jp].x += rtor_0;
              f[jp].y += rtor_1;
              f[jp].z += rtor_2;
            }

            if (EFLAG) {
              evdwl = factor_lj * one_eng;
              sevdwl += evdwl;
              if (eatom) {
                fwtmp += (flt_t)0.5 * evdwl;
                if (NEWTON_PAIR)
                  f[j*2].w += (flt_t)0.5 * evdwl;
              }
            }

            if (NEWTON_PAIR == 0) {
              if (vflag == 1) {
                sv0 += delx_form[jj] * fforce_0;
                sv1 += dely_form[jj] * fforce_1;
                sv2 += delz_form[jj] * fforce_2;
                sv3 += delx_form[jj] * fforce_1;
                sv4 += delx_form[jj] * fforce_2;
                sv5 += dely_form[jj] * fforce_2;
              }
            } // EVFLAG
          #ifdef INTEL_VMASK
          }
          #endif
        } // for jj

        // -------------------------------------------------------------

        if (multiple_forms)
          ierror = 2;

        int ip = i * 2;
        if (NEWTON_PAIR) {
          f[ip].x += fxtmp;
          f[ip].y += fytmp;
          f[ip].z += fztmp;
          ip++;
          f[ip].x += t1tmp;
          f[ip].y += t2tmp;
          f[ip].z += t3tmp;
        } else {
          f[ip].x = fxtmp;
          f[ip].y = fytmp;
          f[ip].z = fztmp;
          ip++;
          f[ip].x = t1tmp;
          f[ip].y = t2tmp;
          f[ip].z = t3tmp;
        }

        if (EFLAG) {
          oevdwl += sevdwl;
          if (eatom) f[i * 2].w += fwtmp;
        }
        if (NEWTON_PAIR == 0) {
          if (vflag == 1) {
            ov0 += sv0;
            ov1 += sv1;
            ov2 += sv2;
            ov3 += sv3;
            ov4 += sv4;
            ov5 += sv5;
          }
        }
      } // for i
      int o_range;
      if (NEWTON_PAIR) {
        o_range = nall;
        if (offload == 0) o_range -= minlocal;
        IP_PRE_omp_range_align(iifrom, iito, tid, o_range, nthreads,
                               sizeof(FORCE_T));
        const int sto = iito * 8;
        const int fst4 = f_stride * 4;
        #if defined(_OPENMP)
        #pragma omp barrier
        #endif
        acc_t *f_scalar = &f_start[0].x;
        acc_t *f_scalar2 = f_scalar + fst4;
        for (int t = 1; t < nthreads; t++) {
          #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #pragma vector aligned
          #endif
          for (int n = iifrom * 8; n < sto; n++)
            f_scalar[n] += f_scalar2[n];
          f_scalar2 += fst4;
        }

        if (vflag==2) {
          const ATOM_T * _noalias const xo = x + minlocal;
          #if defined(LMP_SIMD_COMPILER)
          #pragma novector
          #endif
          for (int n = iifrom; n < iito; n++) {
            const int nt2 = n * 2;
            ov0 += f_start[nt2].x * xo[n].x;
            ov1 += f_start[nt2].y * xo[n].y;
            ov2 += f_start[nt2].z * xo[n].z;
            ov3 += f_start[nt2].y * xo[n].x;
            ov4 += f_start[nt2].z * xo[n].x;
            ov5 += f_start[nt2].z * xo[n].y;
          }
        }
      }

      if (ierror)
        f_start[1].w = ierror;
    } // omp

    if (EFLAG || vflag) {
      if (NEWTON_PAIR == 0) {
        oevdwl *= (acc_t)0.5;
        ov0 *= (acc_t)-0.5;
        ov1 *= (acc_t)-0.5;
        ov2 *= (acc_t)-0.5;
        ov3 *= (acc_t)-0.5;
        ov4 *= (acc_t)-0.5;
        ov5 *= (acc_t)-0.5;
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
  } // offload

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, 2);
  else
    fix->add_result_array(f_start, nullptr, offload, 0, 0, 2);
}

/* ---------------------------------------------------------------------- */

void PairGayBerneIntel::init_style()
{
  PairGayBerne::init_style();
  if (force->newton_pair == 0)
    neighbor->find_request(this)->enable_full();

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  if (force->newton_pair) fix->set_offload_noghost(1);
  _cop = fix->coprocessor_number();
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
void PairGayBerneIntel::pack_force_const(ForceConst<flt_t> &fc,
                                         IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;
  _max_nbors = buffers->get_max_nbors();
  int mthreads = comm->nthreads;
  if (mthreads < buffers->get_off_threads())
    mthreads = buffers->get_off_threads();
  fc.set_ntypes(tp1, _max_nbors, mthreads, memory, _cop);

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

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_lj[0] = 1.0;
  }
  fc.gamma = gamma;
  fc.upsilon = upsilon;
  fc.mu = mu;

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.ijc[i][j].lj1 = lj1[i][j];
      fc.ijc[i][j].lj2 = lj2[i][j];
      fc.ijc[i][j].cutsq = cutsq[i][j];
      fc.ijc[i][j].offset = offset[i][j];
      fc.ijc[i][j].sigma = sigma[i][j];
      fc.ijc[i][j].epsilon = epsilon[i][j];
      fc.ijc[i][j].form = form[i][j];
      fc.ijc[i][j].lshape = lshape[i] * lshape[j];
      fc.lj34[i][j].lj3 = lj3[i][j];
      fc.lj34[i][j].lj4 = lj4[i][j];
    }
    for (int j = 0; j < 4; j++) {
      fc.ic[i].shape2[j] = shape2[i][j];
      fc.ic[i].well[j] = well[i][j];
    }
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  flt_t * special_lj = fc.special_lj;
  FC_PACKED1_T *oijc = fc.ijc[0];
  FC_PACKED2_T *olj34 = fc.lj34[0];
  FC_PACKED3_T *oic = fc.ic;
  int tp1sq = tp1 * tp1;
  if (oijc != nullptr && oic != nullptr) {
    #pragma offload_transfer target(mic:_cop) \
      in(special_lj: length(4) alloc_if(0) free_if(0)) \
      in(oijc,olj34: length(tp1sq) alloc_if(0) free_if(0)) \
      in(oic: length(tp1) alloc_if(0) free_if(0))
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairGayBerneIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                      const int one_length,
                                                      const int nthreads,
                                                      Memory *memory,
                                                      const int cop) {
  if (memory != nullptr) _memory = memory;
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      fc_packed3 *oic = ic;

      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      fc_packed1 *oijc = ijc[0];
      fc_packed2 *olj34 = lj34[0];
      flt_t * orsq_form = rsq_form[0];
      flt_t * odelx_form = delx_form[0];
      flt_t * odely_form = dely_form[0];
      flt_t * odelz_form = delz_form[0];
      int * ojtype_form = jtype_form[0];
      int * ojlist_form = jlist_form[0];

      if (ospecial_lj != nullptr && oijc != nullptr && olj34 != nullptr &&
          orsq_form != nullptr && odelx_form != nullptr && odely_form != nullptr &&
          odelz_form != nullptr && ojtype_form != nullptr && ojlist_form != nullptr &&
          _cop >= 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(ospecial_lj, oijc, olj34, oic: alloc_if(0) free_if(1)) \
          nocopy(orsq_form, odelx_form, odely_form: alloc_if(0) free_if(1)) \
          nocopy(odelz_form, ojtype_form, ojlist_form: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(oic);
      _memory->destroy(ijc);
      _memory->destroy(lj34);
      _memory->destroy(rsq_form);
      _memory->destroy(delx_form);
      _memory->destroy(dely_form);
      _memory->destroy(delz_form);
      _memory->destroy(jtype_form);
      _memory->destroy(jlist_form);
    }

    if (ntypes > 0) {
      _cop = cop;
      _memory->create(ijc, ntypes, ntypes, "fc.ijc");
      _memory->create(lj34, ntypes, ntypes, "fc.lj34");
      _memory->create(ic, ntypes, "fc.ic");
      _memory->create(rsq_form, nthreads, one_length, "rsq_form");
      _memory->create(delx_form, nthreads, one_length, "delx_form");
      _memory->create(dely_form, nthreads, one_length, "dely_form");
      _memory->create(delz_form, nthreads, one_length, "delz_form");
      _memory->create(jtype_form, nthreads, one_length, "jtype_form");
      _memory->create(jlist_form, nthreads, one_length, "jlist_form");

      for (int zn = 0; zn < nthreads; zn++)
        for (int zo = 0; zo < one_length; zo++) {
          rsq_form[zn][zo] = 10.0;
          delx_form[zn][zo] = 10.0;
          dely_form[zn][zo] = 10.0;
          delz_form[zn][zo] = 10.0;
          jtype_form[zn][zo] = 1;
          jlist_form[zn][zo] = 0;
        }

      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      fc_packed1 *oijc = ijc[0];
      fc_packed2 *olj34 = lj34[0];
      fc_packed3 *oic = ic;
      flt_t * orsq_form = rsq_form[0];
      flt_t * odelx_form = delx_form[0];
      flt_t * odely_form = dely_form[0];
      flt_t * odelz_form = delz_form[0];
      int * ojtype_form = jtype_form[0];
      int * ojlist_form = jlist_form[0];
      int off_onel = one_length * nthreads;

      int tp1sq = ntypes*ntypes;
      if (ospecial_lj != nullptr && oijc != nullptr && olj34 != nullptr &&
          oic != nullptr && orsq_form != nullptr && odelx_form != nullptr &&
          odely_form != nullptr && odelz_form != nullptr && ojtype_form !=nullptr &&
          ojlist_form !=nullptr && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj: length(4) alloc_if(1) free_if(0)) \
          nocopy(oijc,olj34: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oic: length(ntypes) alloc_if(1) free_if(0)) \
          in(orsq_form: length(off_onel) alloc_if(1) free_if(0)) \
          in(odelx_form: length(off_onel) alloc_if(1) free_if(0)) \
          in(odely_form: length(off_onel) alloc_if(1) free_if(0)) \
          in(odelz_form: length(off_onel) alloc_if(1) free_if(0)) \
          in(ojtype_form: length(off_onel) alloc_if(1) free_if(0)) \
          in(ojlist_form: length(off_onel) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes = ntypes;
}
