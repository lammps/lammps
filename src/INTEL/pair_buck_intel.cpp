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

#include "pair_buck_intel.h"

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

PairBuckIntel::PairBuckIntel(LAMMPS *lmp) : PairBuck(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

void PairBuckIntel::compute(int eflag, int vflag)
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
void PairBuckIntel::compute(int eflag, int vflag,
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
      eval<1,1>(ovflag, buffers, fc, 0, inum);
    } else {
      eval<1,0>(ovflag, buffers, fc, 0, inum);
    }
  } else {
    if (force->newton_pair) {
      eval<0,1>(ovflag, buffers, fc, 0, inum);
    } else {
      eval<0,0>(ovflag, buffers, fc, 0, inum);
    }
  }
}

/* ---------------------------------------------------------------------- */

template <int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairBuckIntel::eval(const int vflag,
                         IntelBuffers<flt_t,acc_t> *buffers,
                         const ForceConst<flt_t> &fc,
                         const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(nlocal, nall, minlocal);

  ATOM_T * _noalias const x = buffers->get_x();

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;  // NOLINT

  const flt_t * _noalias const special_lj = fc.special_lj;
  const C_FORCE_T * _noalias const c_force = fc.c_force[0];
  const C_ENERGY_T * _noalias const c_energy = fc.c_energy[0];

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int f_stride;
  IP_PRE_get_transfern(NEWTON_PAIR, buffers, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(buffers, fix, tc, f_start, ev_global);

  const int nthreads = tc;
  {


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

      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        const int itype = x[i].w;

        const int ptr_off = itype * ntypes;
        const C_FORCE_T * _noalias const c_forcei = c_force + ptr_off;
        const C_ENERGY_T * _noalias const c_energyi = c_energy + ptr_off;
        const int   * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum);

        acc_t fxtmp,fytmp,fztmp,fwtmp;
        acc_t sevdwl,  sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EFLAG) fwtmp = sevdwl =  (acc_t)0;
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

          flt_t  forcebuck, evdwl;
          forcebuck = evdwl =  (flt_t)0.0;

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
          if (rsq < c_forcei[jtype].cutsq) {
          #endif
            const flt_t r6inv = r2inv * r2inv * r2inv;
            const flt_t rexp = std::exp(-r * c_forcei[jtype].rhoinv);
            forcebuck = r * rexp * c_forcei[jtype].buck1 -
              r6inv * c_forcei[jtype].buck2;

            #ifndef INTEL_VMASK
            if (rsq > c_forcei[jtype].cutsq)
              forcebuck =(flt_t)0.0;
            #endif
            if (EFLAG) {
              evdwl = rexp * c_energyi[jtype].a -
                r6inv * c_energyi[jtype].c -
                c_energyi[jtype].offset;

              #ifndef INTEL_VMASK
              if (rsq > c_forcei[jtype].cutsq)
                evdwl =(flt_t)0.0;
              #endif
            }

            if (sbindex) {
              const flt_t factor_lj = special_lj[sbindex];
              forcebuck *= factor_lj;
              if (EFLAG)
                evdwl *= factor_lj;
            }
            const flt_t fpair =  forcebuck * r2inv;
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
              if (eatom) {
                fwtmp += (flt_t)0.5 * evdwl;
                if (NEWTON_PAIR)
                  f[j].w += (flt_t)0.5 * evdwl;
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
        IP_PRE_ev_tally_atom(NEWTON_PAIR, EFLAG, vflag, f, fwtmp);
      } // for ii

      IP_PRE_fdotr_reduce_omp(NEWTON_PAIR, nall, minlocal, nthreads, f_start,
                              f_stride, x, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
    } // end of omp parallel region

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
      ev_global[1] = (acc_t)0;
      ev_global[2] = ov0;
      ev_global[3] = ov1;
      ev_global[4] = ov2;
      ev_global[5] = ov3;
      ev_global[6] = ov4;
      ev_global[7] = ov5;
    }
  }

  fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, eatom, vflag);
  else
    fix->add_result_array(f_start, nullptr);
}

void PairBuckIntel::init_style()
{
  PairBuck::init_style();

  // augment neighbor list request
  if (force->newton_pair == 0)
    neighbor->find_request(this)->enable_full();

  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  fix->pair_init_check();

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
}

template <class flt_t, class acc_t>
void PairBuckIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;

  fc.set_ntypes(tp1, memory);

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
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.c_force[i][j].buck1 = buck1[i][j];
      fc.c_force[i][j].buck2 = buck2[i][j];
      fc.c_force[i][j].rhoinv = rhoinv[i][j];
      fc.c_force[i][j].cutsq = cutsq[i][j];
      fc.c_energy[i][j].a = a[i][j];
      fc.c_energy[i][j].c = c[i][j];
      fc.c_energy[i][j].offset = offset[i][j];
    }
  }

}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairBuckIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                  Memory *memory) {
  if (memory != nullptr) _memory = memory;
  if ((ntypes != _ntypes )) {
    if (_ntypes > 0) {

      _memory->destroy(c_force);
      _memory->destroy(c_energy);

    }
    if (ntypes > 0) {
      _memory->create(c_force,ntypes,ntypes,"fc.c_force");
      _memory->create(c_energy,ntypes,ntypes,"fc.c_energy");

    }
  }
  _ntypes=ntypes;
}


