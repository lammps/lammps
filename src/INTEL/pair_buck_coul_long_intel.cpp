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
   Contributing author: Rodrigo Canales (RWTH Aachen University)
------------------------------------------------------------------------- */

#include "pair_buck_coul_long_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
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
#define TABLE_T typename ForceConst<flt_t>::table_t

PairBuckCoulLongIntel::PairBuckCoulLongIntel(LAMMPS *lmp) :
  PairBuckCoulLong(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

void PairBuckCoulLongIntel::compute(int eflag, int vflag)
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
void PairBuckCoulLongIntel::compute(int eflag, int vflag,
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

  if (_lrt == 0 && ago != 0 && fix->separate_buffers() == 0) {
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
void PairBuckCoulLongIntel::eval(const int vflag,
                                 IntelBuffers<flt_t,acc_t> *buffers,
                                 const ForceConst<flt_t> &fc,
                                 const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(nlocal, nall, minlocal);

  ATOM_T * _noalias const x = buffers->get_x();
  flt_t * _noalias const q = buffers->get_q();

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;  // NOLINT

  const flt_t * _noalias const special_coul = fc.special_coul;
  const flt_t * _noalias const special_lj = fc.special_lj;
  const flt_t qqrd2e = force->qqrd2e;

  const C_FORCE_T * _noalias const c_force = fc.c_force[0];
  const C_ENERGY_T * _noalias const c_energy = fc.c_energy[0];
  const flt_t * _noalias const rho_inv = fc.rho_inv[0];
  const TABLE_T * _noalias const table = fc.table;
  const flt_t * _noalias const etable = fc.etable;
  const flt_t * _noalias const detable = fc.detable;
  const flt_t * _noalias const ctable = fc.ctable;
  const flt_t * _noalias const dctable = fc.dctable;
  const flt_t g_ewald = fc.g_ewald;
  const flt_t tabinnersq = fc.tabinnersq;

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
  int f_stride;
  IP_PRE_get_transfern(NEWTON_PAIR, buffers, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(buffers, fix, tc, f_start, ev_global);

  const int nthreads = tc;
  {


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
        const C_FORCE_T * _noalias const c_forcei = c_force + ptr_off;
        const C_ENERGY_T * _noalias const c_energyi = c_energy + ptr_off;
        const flt_t * _noalias const rho_invi = rho_inv + ptr_off;

        const int   * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum);

        acc_t fxtmp,fytmp,fztmp,fwtmp;
        acc_t sevdwl, secoul, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const flt_t qtmp = q[i];
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EFLAG) fwtmp = sevdwl = secoul = (acc_t)0;
        if (NEWTON_PAIR == 0)
          if (vflag == 1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

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
          const int jtype = IP_PRE_dword_index(x[j].w);
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < c_forcei[jtype].cutsq) {
            trsq[ej]=rsq;
            tdelx[ej]=delx;
            tdely[ej]=dely;
            tdelz[ej]=delz;
            tjtype[ej]=jtype;
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
          flt_t forcecoul, forcebuck, evdwl, ecoul;
          forcecoul = forcebuck = evdwl = ecoul = (flt_t)0.0;

          const int j = tj[jj] & NEIGHMASK;
          const int sbindex = tj[jj] >> SBBITS & 3;
          const int jtype = tjtype[jj];
          const flt_t rsq = trsq[jj];
          const flt_t r2inv = (flt_t)1.0 / rsq;
          const flt_t r = (flt_t)1.0 / std::sqrt(r2inv);

          #ifdef INTEL_ALLOW_TABLE
          if (!ncoultablebits || rsq <= tabinnersq) {
          #endif
            const flt_t A1 =  0.254829592;
            const flt_t A2 = -0.284496736;
            const flt_t A3 =  1.421413741;
            const flt_t A4 = -1.453152027;
            const flt_t A5 =  1.061405429;
            const flt_t EWALD_F = 1.12837917;
            const flt_t INV_EWALD_P = 1.0 / 0.3275911;

            const flt_t grij = g_ewald * r;
            const flt_t expm2 = std::exp(-grij * grij);
            const flt_t t = INV_EWALD_P / (INV_EWALD_P + grij);
            const flt_t erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            const flt_t prefactor = qqrd2e * qtmp * q[j] / r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
            if (EFLAG) ecoul = prefactor * erfc;

            const flt_t adjust = ((flt_t)1.0 - special_coul[sbindex])*
              prefactor;
            forcecoul -= adjust;
            if (EFLAG) ecoul -= adjust;

          #ifdef INTEL_ALLOW_TABLE
          } else {
            float rsq_lookup = rsq;
            const int itable = (__intel_castf32_u32(rsq_lookup) &
              ncoulmask) >> ncoulshiftbits;
            const flt_t fraction = (rsq_lookup - table[itable].r) *
              table[itable].dr;

            const flt_t tablet = table[itable].f +
              fraction * table[itable].df;
            forcecoul = qtmp * q[j] * tablet;
            if (EFLAG) ecoul = qtmp * q[j] * (etable[itable] +
              fraction * detable[itable]);
            if (sbindex) {
              const flt_t table2 = ctable[itable] +
                fraction * dctable[itable];
              const flt_t prefactor = qtmp * q[j] * table2;
              const flt_t adjust = ((flt_t)1.0 - special_coul[sbindex]) *
                prefactor;
              forcecoul -= adjust;
              if (EFLAG) ecoul -= adjust;
            }
          }
          #endif

          #ifdef INTEL_VMASK
          if (rsq < c_forcei[jtype].cut_ljsq) {
          #endif
            flt_t r6inv = r2inv * r2inv * r2inv;
            flt_t rexp = std::exp(-r * rho_invi[jtype]);
            forcebuck = r * rexp * c_forcei[jtype].buck1 -
              r6inv * c_forcei[jtype].buck2;
            if (EFLAG) evdwl = rexp * c_energyi[jtype].a -
                         r6inv * c_energyi[jtype].c -
                         c_energyi[jtype].offset;

            if (sbindex) {
              const flt_t factor_lj = special_lj[sbindex];
              forcebuck *= factor_lj;
              if (EFLAG) evdwl *= factor_lj;
            }
          #ifdef INTEL_VMASK
          }
          #else
          if (rsq > c_forcei[jtype].cut_ljsq)
            { forcebuck = (flt_t)0.0; evdwl = (flt_t)0.0; }
          #endif

          const flt_t fpair = (forcecoul + forcebuck) * r2inv;
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
            secoul += ecoul;
            if (eatom) {
              fwtmp += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
              if (NEWTON_PAIR)
                f[j].w += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
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
                              f_stride, x, vflag, ov0, ov1, ov2, ov3,
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
  }

  fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, eatom, vflag);
  else
    fix->add_result_array(f_start, nullptr);
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulLongIntel::init_style()
{
  PairBuckCoulLong::init_style();
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

  _lrt = fix->lrt();
}

template <class flt_t, class acc_t>
void PairBuckCoulLongIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int off_ccache = 0;
  buffers->grow_ccache(off_ccache, comm->nthreads, 1);
  _ccache_stride = buffers->ccache_stride();

  int tp1 = atom->ntypes + 1;
  int ntable = 1;
  if (ncoultablebits)
    for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  fc.set_ntypes(tp1, ntable, memory);

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

  fc.g_ewald = force->kspace->g_ewald;
  fc.tabinnersq = tabinnersq;

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_coul[i] = force->special_coul[i];
    fc.special_coul[0] = 1.0;
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      if (cutsq[i][j] < cut_ljsq[i][j])
        error->all(FLERR,
         "Intel variant of lj/buck/coul/long expects lj cutoff<=coulombic");
      fc.c_force[i][j].cutsq = cutsq[i][j];
      fc.c_force[i][j].cut_ljsq = cut_ljsq[i][j];
      fc.c_force[i][j].buck1 = buck1[i][j];
      fc.c_force[i][j].buck2 = buck2[i][j];
      fc.rho_inv[i][j] = rhoinv[i][j];
      fc.c_energy[i][j].a = a[i][j];
      fc.c_energy[i][j].c = c[i][j];
      fc.c_energy[i][j].offset = offset[i][j];
      fc.c_energy[i][j].pad = rhoinv[i][j];
    }
  }

  if (ncoultablebits) {
    for (int i = 0; i < ntable; i++) {
      fc.table[i].r = rtable[i];
      fc.table[i].dr = drtable[i];
      fc.table[i].f = ftable[i];
      fc.table[i].df = dftable[i];
      fc.etable[i] = etable[i];
      fc.detable[i] = detable[i];
      fc.ctable[i] = ctable[i];
      fc.dctable[i] = dctable[i];
    }
  }

}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairBuckCoulLongIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                           const int ntable,
                                                           Memory *memory) {
  if (memory != nullptr) _memory = memory;
  if ((ntypes != _ntypes) || (ntable != _ntable)) {
    if (_ntypes > 0) {

      _memory->destroy(c_force);
      _memory->destroy(c_energy);
      _memory->destroy(table);
      _memory->destroy(rho_inv);
      _memory->destroy(etable);
      _memory->destroy(detable);
      _memory->destroy(ctable);
      _memory->destroy(dctable);
    }
    if (ntypes > 0) {
      _memory->create(c_force,ntypes,ntypes,"fc.c_force");
      _memory->create(c_energy,ntypes,ntypes,"fc.c_energy");
      _memory->create(rho_inv,ntypes,ntypes,"fc.rho_inv");
      _memory->create(table,ntable,"pair:fc.table");
      _memory->create(etable,ntable,"pair:fc.etable");
      _memory->create(detable,ntable,"pair:fc.detable");
      _memory->create(ctable,ntable,"pair:fc.ctable");
      _memory->create(dctable,ntable,"pair:fc.dctable");

    }
  }
  _ntypes=ntypes;
  _ntable=ntable;
}


