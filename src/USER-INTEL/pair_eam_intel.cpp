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
   Contributing authors: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_eam_intel.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "suffix.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2

/* ---------------------------------------------------------------------- */

PairEAMIntel::PairEAMIntel(LAMMPS *lmp) : PairEAM(lmp)
{
  suffix_flag |= Suffix::INTEL;
  fp_float = 0;
}

/* ---------------------------------------------------------------------- */

PairEAMIntel::~PairEAMIntel()
{
  memory->destroy(fp_float);
}

/* ---------------------------------------------------------------------- */

void PairEAMIntel::compute(int eflag, int vflag)
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
void PairEAMIntel::compute(int eflag, int vflag,
                           IntelBuffers<flt_t,acc_t> *buffers,
                           const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  } else evflag = vflag_fdotr = 0;

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
    #pragma omp parallel if(packthreads > 1)
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

/* ---------------------------------------------------------------------- */

template <int ONETYPE, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairEAMIntel::eval(const int offload, const int vflag,
                        IntelBuffers<flt_t,acc_t> *buffers,
                        const ForceConst<flt_t> &fc,
                        const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;

  flt_t *fp_f;
  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    IP_PRE_edge_align(nmax, sizeof(acc_t));
    if (NEWTON_PAIR)
      memory->create(rho,nmax*comm->nthreads,"pair:rho");
    else
      memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    // Use single precision allocation for single/mixed mode
    // Keep double version for single and swap_eam
    if (sizeof(flt_t)==sizeof(float)) {
      memory->destroy(fp_float);
      memory->create(fp_float,nmax,"pair::fp_float");
    }
  }
  if (sizeof(flt_t)==sizeof(float))
    fp_f = (flt_t *)fp_float;
  else
    fp_f = (flt_t *)fp;


  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;
  const FC_PACKED1_T * _noalias const rhor_spline_f = fc.rhor_spline_f;
  const FC_PACKED1_T * _noalias const rhor_spline_e = fc.rhor_spline_e;
  const FC_PACKED2_T * _noalias const z2r_spline_t = fc.z2r_spline_t;
  const FC_PACKED1_T * _noalias const frho_spline_f = fc.frho_spline_f;
  const FC_PACKED1_T * _noalias const frho_spline_e = fc.frho_spline_e;
  const flt_t * _noalias const scale_f = fc.scale_f[0];

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
  int *overflow = fix->get_off_overflow_flag();

  const flt_t frdr = rdr;
  const flt_t frdrho = rdrho;
  const flt_t frhomax = rhomax;
  const flt_t fcutforcesq = cutforcesq;
  const int istride = fc.rhor_istride();
  const int jstride = fc.rhor_jstride();
  const int fstride = fc.frho_stride();

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
      int iifrom, iito, tid;
      IP_PRE_omp_range_id_vec(iifrom, iito, tid, inum, nthreads,
                              INTEL_VECTOR_WIDTH);
      iifrom += astart;
      iito += astart;

      int foff;
      if (NEWTON_PAIR) foff = tid * f_stride - minlocal;
      else foff = -minlocal;
      FORCE_T * _noalias const f = f_start + foff;
      if (NEWTON_PAIR) foff = tid * nmax;
      else foff = 0;
      double * _noalias const trho = rho + foff;
      if (NEWTON_PAIR) {
        memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));
        memset(trho, 0, nall * sizeof(double));
      }

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;

      flt_t oscale;
      int rhor_joff, frho_ioff;
      if (ONETYPE) {
        const int ptr_off=_onetype * ntypes + _onetype;
        oscale = scale_f[ptr_off];
        int rhor_ioff = istride * _onetype;
        rhor_joff = rhor_ioff + _onetype * jstride;
        frho_ioff = fstride * _onetype;
      }
      for (int ii = iifrom; ii < iito; ++ii) {
        const int i = ilist[ii];
        int itype, rhor_ioff;
        if (!ONETYPE) {
          itype = x[i].w;
          rhor_ioff = istride * itype;
        }
        const int * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum, offload);

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;

        acc_t rhoi = (acc_t)0.0;
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
          const flt_t rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < fcutforcesq) {
            trsq[ej]=rsq;
            if (!ONETYPE) tjtype[ej]=x[j].w;
            tj[ej]=jlist[jj];
            ej++;
          }
        }

        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma simd reduction(+:rhoi)
        #endif
        for (int jj = 0; jj < ej; jj++) {
          int jtype;
          const int j = tj[jj] & NEIGHMASK;
          if (!ONETYPE) jtype = tjtype[jj];
          const flt_t rsq = trsq[jj];
          flt_t p = sqrt(rsq)*frdr + (flt_t)1.0;
          int m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,(flt_t)1.0);
          if (!ONETYPE)
            rhor_joff = rhor_ioff + jtype * jstride;
          const int joff = rhor_joff + m;
          flt_t ra;
          ra = ((rhor_spline_e[joff].a*p + rhor_spline_e[joff].b) * p +
                rhor_spline_e[joff].c) * p + rhor_spline_e[joff].d;
          rhoi += ra;
          if (NEWTON_PAIR) {
            if (!ONETYPE) {
              const int ioff = jtype * istride + itype * jstride + m;
              ra = ((rhor_spline_e[ioff].a*p + rhor_spline_e[ioff].b)*p +
                    rhor_spline_e[ioff].c) * p + rhor_spline_e[ioff].d;
            }
            trho[j] += ra;
          }
        } // for jj
        if (NEWTON_PAIR)
          trho[i] += rhoi;
        else
          trho[i] = rhoi;
      } // for i

      #if defined(_OPENMP)
      if (NEWTON_PAIR && nthreads > 1) {
        #pragma omp barrier
        if (tid == 0) {
          const int rcount = nall;
          if (nthreads == 2) {
            double *trho2 = rho + nmax;
            #pragma vector aligned
            #pragma simd
            for (int n = 0; n < rcount; n++)
              rho[n] += trho2[n];
          } else if (nthreads == 4) {
            double *trho2 = rho + nmax;
            double *trho3 = trho2 + nmax;
            double *trho4 = trho3 + nmax;
            #pragma vector aligned
            #pragma simd
            for (int n = 0; n < rcount; n++)
              rho[n] += trho2[n] + trho3[n] + trho4[n];
          } else {
            double *trhon = rho + nmax;
            for (int t = 1; t < nthreads; t++) {
              #pragma vector aligned
              #pragma simd
              for (int n = 0; n < rcount; n++)
                rho[n] += trhon[n];
              trhon += nmax;
            }
          }
        }
      }
      #endif

      // communicate and sum densities

      if (NEWTON_PAIR) {
        if (tid == 0)
          comm->reverse_comm_pair(this);
      }
      #if defined(_OPENMP)
      #pragma omp barrier
      #endif

      // fp = derivative of embedding energy at each atom
      // phi = embedding energy at each atom
      // if rho > rhomax (e.g. due to close approach of two atoms),
      //   will exceed table, so add linear term to conserve energy

      acc_t tevdwl;
      if (EFLAG) tevdwl = (acc_t)0.0;

      #if defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma simd reduction(+:tevdwl)
      #endif
      for (int ii = iifrom; ii < iito; ++ii) {
        const int i = ilist[ii];
        int itype;
        if (!ONETYPE) itype = x[i].w;
        flt_t p = rho[i]*frdrho + (flt_t)1.0;
        int m = static_cast<int> (p);
        m = MAX(1,MIN(m,nrho-1));
        p -= m;
        p = MIN(p,(flt_t)1.0);
        if (!ONETYPE) frho_ioff = itype * fstride;
        const int ioff = frho_ioff + m;
        fp_f[i] = (frho_spline_f[ioff].a*p + frho_spline_f[ioff].b)*p +
          frho_spline_f[ioff].c;
        if (EFLAG) {
          flt_t phi = ((frho_spline_e[ioff].a*p + frho_spline_e[ioff].b)*p +
                       frho_spline_e[ioff].c)*p + frho_spline_e[ioff].d;
          if (rho[i] > frhomax) phi += fp_f[i] * (rho[i]-frhomax);
          if (!ONETYPE) {
            const int ptr_off=itype*ntypes + itype;
            oscale = scale_f[ptr_off];
          }
          phi *= oscale;
          tevdwl += phi;
          if (eatom) f[i].w += phi;
        }
      }
      if (EFLAG) oevdwl += tevdwl;


      // communicate derivative of embedding function

      #if defined(_OPENMP)
      #pragma omp barrier
      #endif

      if (tid == 0)
        comm->forward_comm_pair(this);
      if (NEWTON_PAIR) memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      #if defined(_OPENMP)
      #pragma omp barrier
      #endif

      // compute forces on each atom
      // loop over neighbors of my atoms

      for (int ii = iifrom; ii < iito; ++ii) {
        const int i = ilist[ii];
        int itype, rhor_ioff;
        const flt_t * _noalias scale_fi;
        if (!ONETYPE) {
          itype = x[i].w;
          rhor_ioff = istride * itype;
          scale_fi = scale_f + itype*ntypes;
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
          if (vflag==1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

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
          const flt_t rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < fcutforcesq) {
            trsq[ej]=rsq;
            tdelx[ej]=delx;
            tdely[ej]=dely;
            tdelz[ej]=delz;
            if (!ONETYPE) tjtype[ej]=x[j].w;
            tj[ej]=jlist[jj];
            ej++;
          }
        }

        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                                 sv0, sv1, sv2, sv3, sv4, sv5)
        #endif
        for (int jj = 0; jj < ej; jj++) {
          int jtype;
          const int j = tj[jj] & NEIGHMASK;
          if (!ONETYPE) jtype = tjtype[jj];
          const flt_t rsq = trsq[jj];
          const flt_t r = sqrt(rsq);
          flt_t p = r*frdr + (flt_t)1.0;
          int m = static_cast<int> (p);
          m = MIN(m,nr-1);
          p -= m;
          p = MIN(p,(flt_t)1.0);
          if (!ONETYPE)
            rhor_joff = rhor_ioff + jtype * jstride;
          const int joff = rhor_joff + m;
          const flt_t rhojp = (rhor_spline_f[joff].a*p +
                               rhor_spline_f[joff].b)*p +
            rhor_spline_f[joff].c;
          flt_t rhoip;
          if (!ONETYPE) {
            const int ioff = jtype * istride + itype * jstride + m;
            rhoip = (rhor_spline_f[ioff].a*p + rhor_spline_f[ioff].b)*p +
              rhor_spline_f[ioff].c;
          } else
            rhoip = rhojp;
          const flt_t z2p = (z2r_spline_t[joff].a*p +
                             z2r_spline_t[joff].b)*p +
            z2r_spline_t[joff].c;
          const flt_t z2 = ((z2r_spline_t[joff].d*p +
                             z2r_spline_t[joff].e)*p +
                            z2r_spline_t[joff].f)*p +
            z2r_spline_t[joff].g;

          const flt_t recip = (flt_t)1.0/r;
          const flt_t phi = z2*recip;
          const flt_t phip = z2p*recip - phi*recip;
          const flt_t psip = fp_f[i]*rhojp + fp_f[j]*rhoip + phip;
          if (!ONETYPE)
            oscale = scale_fi[jtype];
          const flt_t fpair = -oscale*psip*recip;

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
            const flt_t evdwl = oscale*phi;
            sevdwl += evdwl;
            if (eatom) {
              fwtmp += (flt_t)0.5 * evdwl;
              if (NEWTON_PAIR)
                f[j].w += (flt_t)0.5 * evdwl;
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
          sevdwl *= (acc_t)0.5;
        }

        IP_PRE_ev_tally_atom(NEWTON_PAIR, EFLAG, vflag, f, fwtmp);
      } // for i

      IP_PRE_fdotr_reduce_omp(NEWTON_PAIR, nall, minlocal, nthreads, f_start,
                              f_stride, x, offload, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
    } /// omp

    IP_PRE_fdotr_reduce(NEWTON_PAIR, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      if (NEWTON_PAIR == 0) {
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
  }

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMIntel::init_style()
{
  PairEAM::init_style();
  if (force->newton_pair == 0) {
    neighbor->requests[neighbor->nrequest-1]->half = 0;
    neighbor->requests[neighbor->nrequest-1]->full = 1;
  }
  neighbor->requests[neighbor->nrequest-1]->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_balance() != 0.0)
    error->all(FLERR,
      "Offload for eam/intel is not yet available. Set balance to 0.");
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
void PairEAMIntel::pack_force_const(ForceConst<flt_t> &fc,
                                    IntelBuffers<flt_t,acc_t> *buffers)
{
  int off_ccache = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop >= 0) off_ccache = 1;
  #endif
  buffers->grow_ccache(off_ccache, comm->nthreads, 1);
  _ccache_stride = buffers->ccache_stride();

  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,nr,nrho,memory,_cop);

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

  _onetype=-1;
  double oldscale=-1;
  for (int i = 1; i < tp1; i++) {
    int ioff = i * fc.frho_stride();
    for (int k = 0; k < nrho + 1; k++) {
      fc.frho_spline_f[ioff + k].a = frho_spline[type2frho[i]][k][0];
      fc.frho_spline_f[ioff + k].b = frho_spline[type2frho[i]][k][1];
      fc.frho_spline_f[ioff + k].c = frho_spline[type2frho[i]][k][2];
      fc.frho_spline_e[ioff + k].a = frho_spline[type2frho[i]][k][3];
      fc.frho_spline_e[ioff + k].b = frho_spline[type2frho[i]][k][4];
      fc.frho_spline_e[ioff + k].c = frho_spline[type2frho[i]][k][5];
      fc.frho_spline_e[ioff + k].d = frho_spline[type2frho[i]][k][6];
    }
    ioff = i * fc.rhor_istride();
    for (int j = 1; j < tp1; j++) {
      fc.scale_f[i][j] = scale[i][j];
      if (type2rhor[i][j] >= 0) {
        const int joff = ioff + j * fc.rhor_jstride();
        for (int k = 0; k < nr + 1; k++) {
          if (type2rhor[j][i] != type2rhor[i][j])
            _onetype = 0;
          else if (_onetype < 0)
            _onetype = i;
          if (oldscale < 0)
            oldscale = scale[i][j];
          else
            if (oldscale != scale[i][j])
              _onetype = 0;
          fc.rhor_spline_f[joff + k].a=rhor_spline[type2rhor[j][i]][k][0];
          fc.rhor_spline_f[joff + k].b=rhor_spline[type2rhor[j][i]][k][1];
          fc.rhor_spline_f[joff + k].c=rhor_spline[type2rhor[j][i]][k][2];
          fc.rhor_spline_e[joff + k].a=rhor_spline[type2rhor[j][i]][k][3];
          fc.rhor_spline_e[joff + k].b=rhor_spline[type2rhor[j][i]][k][4];
          fc.rhor_spline_e[joff + k].c=rhor_spline[type2rhor[j][i]][k][5];
          fc.rhor_spline_e[joff + k].d=rhor_spline[type2rhor[j][i]][k][6];
          fc.z2r_spline_t[joff + k].a=z2r_spline[type2z2r[j][i]][k][0];
          fc.z2r_spline_t[joff + k].b=z2r_spline[type2z2r[j][i]][k][1];
          fc.z2r_spline_t[joff + k].c=z2r_spline[type2z2r[j][i]][k][2];
          fc.z2r_spline_t[joff + k].d=z2r_spline[type2z2r[j][i]][k][3];
          fc.z2r_spline_t[joff + k].e=z2r_spline[type2z2r[j][i]][k][4];
          fc.z2r_spline_t[joff + k].f=z2r_spline[type2z2r[j][i]][k][5];
          fc.z2r_spline_t[joff + k].g=z2r_spline[type2z2r[j][i]][k][6];
        }
      }
    }
  }
  if (_onetype < 0) _onetype = 0;
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairEAMIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                 const int nr, const int nrho,
                                                 Memory *memory,
                                                 const int cop) {
  if (ntypes != _ntypes || nr + 1 > _nr || nrho + 1 > _nrho) {
    if (_ntypes > 0) {
      _memory->destroy(rhor_spline_f);
      _memory->destroy(rhor_spline_e);
      _memory->destroy(frho_spline_f);
      _memory->destroy(frho_spline_e);
      _memory->destroy(z2r_spline_t);
      _memory->destroy(scale_f);
    }
    if (ntypes > 0) {
      _cop = cop;
      _nr = nr + 1;
      IP_PRE_edge_align(_nr, sizeof(flt_t));
      memory->create(rhor_spline_f,ntypes*ntypes*_nr,"fc.rhor_spline_f");
      memory->create(rhor_spline_e,ntypes*ntypes*_nr,"fc.rhor_spline_e");
      memory->create(z2r_spline_t,ntypes*ntypes*_nr,"fc.z2r_spline_t");
      _nrho = nrho + 1;
      IP_PRE_edge_align(_nrho, sizeof(flt_t));
      memory->create(frho_spline_f,ntypes*_nrho,"fc.frho_spline_f");
      memory->create(frho_spline_e,ntypes*_nrho,"fc.frho_spline_e");
      memory->create(scale_f,ntypes,ntypes,"fc.scale_f");
    }
  }
  _ntypes = ntypes;
  _memory = memory;
}

/* ---------------------------------------------------------------------- */

int PairEAMIntel::pack_forward_comm(int n, int *list, double *buf,
                                    int /*pbc_flag*/, int * /*pbc*/)
{
  if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    return pack_forward_comm(n, list, buf, fp);
  else
    return pack_forward_comm(n, list, buf, fp_float);
}

/* ---------------------------------------------------------------------- */

void PairEAMIntel::unpack_forward_comm(int n, int first, double *buf)
{
  if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    unpack_forward_comm(n, first, buf, fp);
  else
    unpack_forward_comm(n, first, buf, fp_float);
}

/* ---------------------------------------------------------------------- */

template<class flt_t>
int PairEAMIntel::pack_forward_comm(int n, int *list, double *buf,
                                    flt_t *fp_f)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp_f[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class flt_t>
void PairEAMIntel::unpack_forward_comm(int n, int first, double *buf,
                                       flt_t *fp_f)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp_f[i] = buf[m++];
}

