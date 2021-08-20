// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "omp_compat.h"

#include <cmath>
#include "dihedral_charmm_intel.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "error.h"

#ifdef LMP_USE_AVXCD
#if (__INTEL_COMPILER_BUILD_DATE > 20160414)
#define LMP_USE_AVXCD_DHC
#endif
#endif

#ifdef LMP_USE_AVXCD_DHC
#include "intel_simd.h"
using namespace ip_simd;
#endif

#include "suffix.h"
using namespace LAMMPS_NS;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

DihedralCharmmIntel::DihedralCharmmIntel(class LAMMPS *lmp)
  : DihedralCharmm(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void DihedralCharmmIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    DihedralCharmm::compute(eflag, vflag);
    return;
  }
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void DihedralCharmmIntel::compute(int eflag, int vflag,
                                  IntelBuffers<flt_t,acc_t> *buffers,
                                  const ForceConst<flt_t> &fc)
{
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"INTEL package does not support per-atom stress");

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (weightflag && vflag_global == VIRIAL_FDOTR)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  if (evflag) {
    if (vflag && !eflag) {
      if (force->newton_bond)
        eval<0,1,1>(vflag, buffers, fc);
      else
        eval<0,1,0>(vflag, buffers, fc);
    } else {
      if (force->newton_bond)
        eval<1,1,1>(vflag, buffers, fc);
      else
        eval<1,1,0>(vflag, buffers, fc);
    }
  } else {
    if (force->newton_bond)
      eval<0,0,1>(vflag, buffers, fc);
    else
      eval<0,0,0>(vflag, buffers, fc);
  }
}

#ifndef LMP_USE_AVXCD_DHC

template <int EFLAG, int VFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void DihedralCharmmIntel::eval(const int vflag,
                               IntelBuffers<flt_t,acc_t> *buffers,
                               const ForceConst<flt_t> &fc)

{
  const int inum = neighbor->ndihedrallist;
  if (inum == 0) return;

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  int f_stride;
  if (NEWTON_BOND) f_stride = buffers->get_stride(nall);
  else f_stride = buffers->get_stride(nlocal);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(0, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  acc_t oedihedral, ov0, ov1, ov2, ov3, ov4, ov5;
  acc_t oevdwl, oecoul, opv0, opv1, opv2, opv3, opv4, opv5;
  if (EFLAG) oevdwl = oecoul = oedihedral = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
    opv0 = opv1 = opv2 = opv3 = opv4 = opv5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(f_start,f_stride,fc)           \
    reduction(+:oevdwl,oecoul,oedihedral,ov0,ov1,ov2,ov3,ov4,ov5, \
              opv0,opv1,opv2,opv3,opv4,opv5)
  #endif
  {
    #if defined(LMP_SIMD_COMPILER_TEST)
    int nfrom, nto, tid;
    IP_PRE_omp_range_id(nfrom, nto, tid, inum, nthreads);
    #else
    int nfrom, npl, nto, tid;
    IP_PRE_omp_stride_id(nfrom, npl, nto, tid, inum, nthreads);
    #endif

    FORCE_T * _noalias const f = f_start + (tid * f_stride);
    if (fix->need_zero(tid))
      memset(f, 0, f_stride * sizeof(FORCE_T));

    const int5_t * _noalias const dihedrallist =
      (int5_t *) neighbor->dihedrallist[0];
    const flt_t qqrd2e = force->qqrd2e;

    acc_t sedihedral, sv0, sv1, sv2, sv3, sv4, sv5;
    acc_t sevdwl, secoul, spv0, spv1, spv2, spv3, spv4, spv5;
    if (EFLAG) sevdwl = secoul = sedihedral = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
      spv0 = spv1 = spv2 = spv3 = spv4 = spv5 = (acc_t)0.0;
    }

    #if defined(LMP_SIMD_COMPILER_TEST)
#if defined(USE_OMP_SIMD)
    #pragma omp simd reduction(+:sedihedral, sevdwl, secoul, sv0, sv1, sv2, \
                               sv3, sv4, sv5, spv0, spv1, spv2, spv3, spv4, \
                               spv5)
#else
    #pragma simd reduction(+:sedihedral, sevdwl, secoul, sv0, sv1, sv2, \
                           sv3, sv4, sv5, spv0, spv1, spv2, spv3, spv4, \
                           spv5)
#endif
    #pragma vector aligned
    for (int n = nfrom; n < nto; n++) {
    #endif
    for (int n = nfrom; n < nto; n += npl) {
      const int i1 = dihedrallist[n].a;
      const int i2 = dihedrallist[n].b;
      const int i3 = dihedrallist[n].c;
      const int i4 = dihedrallist[n].d;
      const int type = dihedrallist[n].t;

      // 1st bond

      const flt_t vb1x = x[i1].x - x[i2].x;
      const flt_t vb1y = x[i1].y - x[i2].y;
      const flt_t vb1z = x[i1].z - x[i2].z;
      const int itype = x[i1].w;

      // 2nd bond

      const flt_t vb2xm = x[i2].x - x[i3].x;
      const flt_t vb2ym = x[i2].y - x[i3].y;
      const flt_t vb2zm = x[i2].z - x[i3].z;

      // 3rd bond

      const flt_t vb3x = x[i4].x - x[i3].x;
      const flt_t vb3y = x[i4].y - x[i3].y;
      const flt_t vb3z = x[i4].z - x[i3].z;
      const int jtype = x[i4].w;

      // 1-4

      const flt_t delx = x[i1].x - x[i4].x;
      const flt_t dely = x[i1].y - x[i4].y;
      const flt_t delz = x[i1].z - x[i4].z;


      // c,s calculation

      const flt_t ax = vb1y*vb2zm - vb1z*vb2ym;
      const flt_t ay = vb1z*vb2xm - vb1x*vb2zm;
      const flt_t az = vb1x*vb2ym - vb1y*vb2xm;
      const flt_t bx = vb3y*vb2zm - vb3z*vb2ym;
      const flt_t by = vb3z*vb2xm - vb3x*vb2zm;
      const flt_t bz = vb3x*vb2ym - vb3y*vb2xm;

      const flt_t rasq = ax*ax + ay*ay + az*az;
      const flt_t rbsq = bx*bx + by*by + bz*bz;
      const flt_t rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      const flt_t rg = sqrt(rgsq);

      flt_t rginv, ra2inv, rb2inv;
      rginv = ra2inv = rb2inv = (flt_t)0.0;
      if (rg > 0) rginv = (flt_t)1.0/rg;
      if (rasq > 0) ra2inv = (flt_t)1.0/rasq;
      if (rbsq > 0) rb2inv = (flt_t)1.0/rbsq;
      const flt_t rabinv = sqrt(ra2inv*rb2inv);

      flt_t c = (ax*bx + ay*by + az*bz)*rabinv;
      const flt_t s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      // error check
      #ifndef LMP_SIMD_COMPILER_TEST
      if (c > PTOLERANCE || c < MTOLERANCE)
        problem(FLERR, i1, i2, i3, i4);
      #endif

      if (c > (flt_t)1.0) c = (flt_t)1.0;
      if (c < (flt_t)-1.0) c = (flt_t)-1.0;

      const flt_t tcos_shift = fc.bp[type].cos_shift;
      const flt_t tsin_shift = fc.bp[type].sin_shift;
      const flt_t tk = fc.bp[type].k;
      const int m = fc.bp[type].multiplicity;

      flt_t p = (flt_t)1.0;
      flt_t ddf1, df1;
      ddf1 = df1 = (flt_t)0.0;

      for (int i = 0; i < m; i++) {
        ddf1 = p*c - df1*s;
        df1 = p*s + df1*c;
        p = ddf1;
      }

      p = p*tcos_shift + df1*tsin_shift;
      df1 = df1*tcos_shift - ddf1*tsin_shift;
      df1 *= -m;
      p += (flt_t)1.0;

      if (m == 0) {
        p = (flt_t)1.0 + tcos_shift;
        df1 = (flt_t)0.0;
      }

      const flt_t fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
      const flt_t hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
      const flt_t fga = fg*ra2inv*rginv;
      const flt_t hgb = hg*rb2inv*rginv;
      const flt_t gaa = -ra2inv*rg;
      const flt_t gbb = rb2inv*rg;

      const flt_t dtfx = gaa*ax;
      const flt_t dtfy = gaa*ay;
      const flt_t dtfz = gaa*az;
      const flt_t dtgx = fga*ax - hgb*bx;
      const flt_t dtgy = fga*ay - hgb*by;
      const flt_t dtgz = fga*az - hgb*bz;
      const flt_t dthx = gbb*bx;
      const flt_t dthy = gbb*by;
      const flt_t dthz = gbb*bz;

      const flt_t df = -tk * df1;

      const flt_t sx2 = df*dtgx;
      const flt_t sy2 = df*dtgy;
      const flt_t sz2 = df*dtgz;

      flt_t f1x = df*dtfx;
      flt_t f1y = df*dtfy;
      flt_t f1z = df*dtfz;

      const flt_t f2x = sx2 - f1x;
      const flt_t f2y = sy2 - f1y;
      const flt_t f2z = sz2 - f1z;

      flt_t f4x = df*dthx;
      flt_t f4y = df*dthy;
      flt_t f4z = df*dthz;

      const flt_t f3x = -sx2 - f4x;
      const flt_t f3y = -sy2 - f4y;
      const flt_t f3z = -sz2 - f4z;

      if (EFLAG || VFLAG) {
        flt_t deng;
        if (EFLAG) deng = tk * p;
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, deng, i1, i2, i3,
                              i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y,
                              f4z, vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm,
                              vb3x, vb3y, vb3z, sedihedral, f, NEWTON_BOND,
                              nlocal, sv0, sv1, sv2, sv3, sv4, sv5);
      }


      #if defined(LMP_SIMD_COMPILER_TEST)
#if defined(USE_OMP_SIMD)
      #pragma omp ordered simd
#else
      #pragma simdoff
#endif
      #endif
      {
        if (NEWTON_BOND || i2 < nlocal) {
          f[i2].x += f2x;
          f[i2].y += f2y;
          f[i2].z += f2z;
        }

        if (NEWTON_BOND || i3 < nlocal) {
          f[i3].x += f3x;
          f[i3].y += f3y;
          f[i3].z += f3z;
        }
      }

      // 1-4 LJ and Coulomb interactions
      // tally energy/virial in pair, using newton_bond as newton flag

      const flt_t tweight = fc.weight[type];
      const flt_t rsq = delx*delx + dely*dely + delz*delz;
      const flt_t r2inv = (flt_t)1.0/rsq;
      const flt_t r6inv = r2inv*r2inv*r2inv;

      flt_t forcecoul;
      if (implicit) forcecoul = qqrd2e * q[i1]*q[i4]*r2inv;
      else forcecoul = qqrd2e * q[i1]*q[i4]*sqrt(r2inv);
      const flt_t forcelj = r6inv * (fc.ljp[itype][jtype].lj1*r6inv -
                                     fc.ljp[itype][jtype].lj2);
      const flt_t fpair = tweight * (forcelj+forcecoul)*r2inv;

      if (NEWTON_BOND || i1 < nlocal) {
        f1x += delx*fpair;
        f1y += dely*fpair;
        f1z += delz*fpair;
      }
      if (NEWTON_BOND || i4 < nlocal) {
        f4x -= delx*fpair;
        f4y -= dely*fpair;
        f4z -= delz*fpair;
      }

      if (EFLAG || VFLAG) {
        flt_t ev_pre = (flt_t)0;
        if (NEWTON_BOND || i1 < nlocal)
          ev_pre += (flt_t)0.5;
        if (NEWTON_BOND || i4 < nlocal)
          ev_pre += (flt_t)0.5;

        if (EFLAG) {
          flt_t ecoul, evdwl;
          ecoul = tweight * forcecoul;
          evdwl = tweight * r6inv * (fc.ljp[itype][jtype].lj3*r6inv -
                                     fc.ljp[itype][jtype].lj4);
          secoul += ev_pre * ecoul;
          sevdwl += ev_pre * evdwl;
          if (eatom) {
            evdwl *= (flt_t)0.5;
            evdwl += (flt_t)0.5 * ecoul;
            if (NEWTON_BOND || i1 < nlocal)
              f[i1].w += evdwl;
            if (NEWTON_BOND || i4 < nlocal)
              f[i4].w += evdwl;
          }
        }
        //            IP_PRE_ev_tally_nbor(vflag, ev_pre, fpair,
        //                                 delx, dely, delz);
        if (VFLAG && vflag) {
          spv0 += ev_pre * delx * delx * fpair;
          spv1 += ev_pre * dely * dely * fpair;
          spv2 += ev_pre * delz * delz * fpair;
          spv3 += ev_pre * delx * dely * fpair;
          spv4 += ev_pre * delx * delz * fpair;
          spv5 += ev_pre * dely * delz * fpair;
        }
      }

      // apply force to each of 4 atoms
      #if defined(LMP_SIMD_COMPILER_TEST)
#if defined(USE_OMP_SIMD)
      #pragma omp ordered simd
#else
      #pragma simdoff
#endif
      #endif
      {
        if (NEWTON_BOND || i1 < nlocal) {
          f[i1].x += f1x;
          f[i1].y += f1y;
          f[i1].z += f1z;
        }

        if (NEWTON_BOND || i4 < nlocal) {
          f[i4].x += f4x;
          f[i4].y += f4y;
          f[i4].z += f4z;
        }
      }
    } // for n
    if (EFLAG) {
      oedihedral += sedihedral;
      oecoul += secoul;
      oevdwl += sevdwl;
    }
    if (VFLAG && vflag) {
      ov0 += sv0; ov1 += sv1; ov2 += sv2; ov3 += sv3; ov4 += sv4; ov5 += sv5;
      opv0 += spv0; opv1 += spv1; opv2 += spv2;
      opv3 += spv3; opv4 += spv4; opv5 += spv5;
    }
  } // omp parallel

  if (EFLAG) {
    energy += oedihedral;
    force->pair->eng_vdwl += oevdwl;
    force->pair->eng_coul += oecoul;
  }
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
    force->pair->virial[0] += opv0;
    force->pair->virial[1] += opv1;
    force->pair->virial[2] += opv2;
    force->pair->virial[3] += opv3;
    force->pair->virial[4] += opv4;
    force->pair->virial[5] += opv5;
  }

  fix->set_reduce_flag();
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

template <int EFLAG, int VFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void DihedralCharmmIntel::eval(const int vflag,
                               IntelBuffers<flt_t,acc_t> *buffers,
                               const ForceConst<flt_t> &fc)

{
  typedef typename SIMD_type<flt_t>::SIMD_vec SIMD_flt_t;
  typedef typename SIMD_type<acc_t>::SIMD_vec SIMD_acc_t;
  const int swidth = SIMD_type<flt_t>::width();

  const int inum = neighbor->ndihedrallist;
  if (inum == 0) return;

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  int f_stride;
  if (NEWTON_BOND) f_stride = buffers->get_stride(nall);
  else f_stride = buffers->get_stride(nlocal);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(0, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  acc_t oedihedral, ov0, ov1, ov2, ov3, ov4, ov5;
  acc_t oevdwl, oecoul, opv0, opv1, opv2, opv3, opv4, opv5;
  if (EFLAG) oevdwl = oecoul = oedihedral = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
    opv0 = opv1 = opv2 = opv3 = opv4 = opv5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(f_start,f_stride,fc)           \
    reduction(+:oevdwl,oecoul,oedihedral,ov0,ov1,ov2,ov3,ov4,ov5, \
              opv0,opv1,opv2,opv3,opv4,opv5)
  #endif
  {
    int nfrom, npl, nto, tid;
    IP_PRE_omp_stride_id_vec(nfrom, npl, nto, tid, inum, nthreads,
                             swidth);

    FORCE_T * _noalias const f = f_start + (tid * f_stride);
    if (fix->need_zero(tid))
      memset(f, 0, f_stride * sizeof(FORCE_T));

    const int * _noalias const dihedrallist =
      (int *) neighbor->dihedrallist[0];
    const flt_t * _noalias const weight = &(fc.weight[0]);
    const flt_t * _noalias const x_f = &(x[0].x);
    const flt_t * _noalias const cos_shift = &(fc.bp[0].cos_shift);
    const flt_t * _noalias const sin_shift = &(fc.bp[0].sin_shift);
    const flt_t * _noalias const k = &(fc.bp[0].k);
    const int * _noalias const multiplicity = &(fc.bp[0].multiplicity);
    const flt_t * _noalias const plj1 = &(fc.ljp[0][0].lj1);
    const flt_t * _noalias const plj2 = &(fc.ljp[0][0].lj2);
    const flt_t * _noalias const plj3 = &(fc.ljp[0][0].lj3);
    const flt_t * _noalias const plj4 = &(fc.ljp[0][0].lj4);
    acc_t * _noalias const pforce= &(f[0].x);
    acc_t * _noalias const featom = &(f[0].w);
    const flt_t qqrd2e = force->qqrd2e;

    SIMD_acc_t sedihedral, sv0, sv1, sv2, sv3, sv4, sv5;
    SIMD_acc_t sevdwl, secoul, spv0, spv1, spv2, spv3, spv4, spv5;
    if (EFLAG) {
      sevdwl = SIMD_set((acc_t)0.0);
      secoul = SIMD_set((acc_t)0.0);
      sedihedral = SIMD_set((acc_t)0.0);
    }
    if (VFLAG && vflag) {
      sv0 = SIMD_set((acc_t)0.0);
      sv1 = SIMD_set((acc_t)0.0);
      sv2 = SIMD_set((acc_t)0.0);
      sv3 = SIMD_set((acc_t)0.0);
      sv4 = SIMD_set((acc_t)0.0);
      sv5 = SIMD_set((acc_t)0.0);
      spv0 = SIMD_set((acc_t)0.0);
      spv1 = SIMD_set((acc_t)0.0);
      spv2 = SIMD_set((acc_t)0.0);
      spv3 = SIMD_set((acc_t)0.0);
      spv4 = SIMD_set((acc_t)0.0);
      spv5 = SIMD_set((acc_t)0.0);
    }

    SIMD_int n_offset = SIMD_set(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                                 55, 60, 65, 70, 75) + (nfrom * 5);
    const int nto5 = nto * 5;
    const int nlocals4 = nlocal << 4;
    const SIMD_int simd_nlocals4 = SIMD_set(nlocals4);
    const int ntypes = atom->ntypes + 1;

    for (int n = nfrom; n < nto; n += npl) {
      SIMD_mask nmask = n_offset < nto5;
      SIMD_int i1 = SIMD_gather(nmask, dihedrallist, n_offset);
      const SIMD_flt_t q1 = SIMD_gather(nmask, q, i1);
      i1 = i1 << 4;
      const SIMD_int i2 = SIMD_gather(nmask, dihedrallist+1, n_offset) << 4;
      const SIMD_int i3 = SIMD_gather(nmask, dihedrallist+2, n_offset) << 4;
      SIMD_int i4 = SIMD_gather(nmask, dihedrallist+3, n_offset);
      const SIMD_flt_t q4 = SIMD_gather(nmask, q, i4);
      i4 = i4 << 4;
      SIMD_int type = SIMD_gather(nmask, dihedrallist+4, n_offset);
      const SIMD_flt_t tweight = SIMD_gather(nmask, weight, type);
      type = type << 2;
      n_offset = n_offset + npl * 5;

      // 1st bond

      SIMD_flt_t x1, x2, y1, y2, z1, z2;
      SIMD_int itype;

      SIMD_atom_gather(nmask, x_f, i1, x1, y1, z1, itype);
      SIMD_atom_gather(nmask, x_f, i2, x2, y2, z2);

      const SIMD_flt_t vb1x = x1 - x2;
      const SIMD_flt_t vb1y = y1 - y2;
      const SIMD_flt_t vb1z = z1 - z2;

      // 2nd bond

      SIMD_flt_t x3, y3, z3;

      SIMD_atom_gather(nmask, x_f, i3, x3, y3, z3);

      const SIMD_flt_t vb2xm = x2 - x3;
      const SIMD_flt_t vb2ym = y2 - y3;
      const SIMD_flt_t vb2zm = z2 - z3;

      // 3rd bond

      SIMD_flt_t x4, y4, z4;
      SIMD_int jtype;

      SIMD_atom_gather(nmask, x_f, i4, x4, y4, z4, jtype);

      const SIMD_flt_t vb3x = x4 - x3;
      const SIMD_flt_t vb3y = y4 - y3;
      const SIMD_flt_t vb3z = z4 - z3;

      // 1-4

      const SIMD_flt_t delx = x1 - x4;
      const SIMD_flt_t dely = y1 - y4;
      const SIMD_flt_t delz = z1 - z4;

      // c,s calculation

      const SIMD_flt_t ax = vb1y*vb2zm - vb1z*vb2ym;
      const SIMD_flt_t ay = vb1z*vb2xm - vb1x*vb2zm;
      const SIMD_flt_t az = vb1x*vb2ym - vb1y*vb2xm;
      const SIMD_flt_t bx = vb3y*vb2zm - vb3z*vb2ym;
      const SIMD_flt_t by = vb3z*vb2xm - vb3x*vb2zm;
      const SIMD_flt_t bz = vb3x*vb2ym - vb3y*vb2xm;

      const SIMD_flt_t rasq = ax*ax + ay*ay + az*az;
      const SIMD_flt_t rbsq = bx*bx + by*by + bz*bz;
      const SIMD_flt_t rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      const SIMD_flt_t rg = SIMD_sqrt(rgsq);

      const SIMD_flt_t szero = SIMD_set((flt_t)0.0);
      const SIMD_flt_t rginv = SIMD_rcpz(rg > szero, rg);
      const SIMD_flt_t ra2inv = SIMD_rcpz(rasq > szero, rasq);
      const SIMD_flt_t rb2inv = SIMD_rcpz(rbsq > szero, rbsq);
      const SIMD_flt_t rabinv = SIMD_sqrt(ra2inv*rb2inv);

      SIMD_flt_t c = (ax*bx + ay*by + az*bz)*rabinv;
      const SIMD_flt_t s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      // error check
      const SIMD_flt_t one = SIMD_set((flt_t)1.0);
      const SIMD_flt_t mone = SIMD_set((flt_t)-1.0);

      const SIMD_flt_t ptol = SIMD_set(PTOLERANCE);
      const SIMD_flt_t ntol = SIMD_set(MTOLERANCE);
      if (c > ptol || c < ntol)
        if (screen)
          error->warning(FLERR,"Dihedral problem.");

      c = SIMD_set(c, c > one, one);
      c = SIMD_set(c, c < mone, mone);

      const SIMD_flt_t tcos_shift = SIMD_gather(nmask, cos_shift, type);
      const SIMD_flt_t tsin_shift = SIMD_gather(nmask, sin_shift, type);
      const SIMD_flt_t tk = SIMD_gather(nmask, k, type);
      const SIMD_int m = SIMD_gatherz_offset<flt_t>(nmask, multiplicity, type);

      SIMD_flt_t p(one);
      SIMD_flt_t ddf1(szero);
      SIMD_flt_t df1(szero);

      const int m_max = SIMD_max(m);

      for (int i = 0; i < m_max; i++) {
        const SIMD_mask my_m = i < m;
        ddf1 = SIMD_set(ddf1, my_m, p*c - df1*s);
        df1 = SIMD_set(df1, my_m, p*s + df1*c);
        p = SIMD_set(p, my_m, ddf1);
      }

      SIMD_flt_t multf;
      SIMD_cast(-m,multf);
      p = p*tcos_shift + df1*tsin_shift;
      df1 = df1*tcos_shift - ddf1*tsin_shift;
      df1 = df1 * multf;
      p = p + one;

      SIMD_mask mzero = (m == SIMD_set((int)0));
      p = SIMD_set(p, mzero, one + tcos_shift);
      df1 = SIMD_set(df1, mzero, szero);

      const SIMD_flt_t fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
      const SIMD_flt_t hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
      const SIMD_flt_t fga = fg*ra2inv*rginv;
      const SIMD_flt_t hgb = hg*rb2inv*rginv;
      const SIMD_flt_t gaa = -ra2inv*rg;
      const SIMD_flt_t gbb = rb2inv*rg;

      const SIMD_flt_t dtfx = gaa*ax;
      const SIMD_flt_t dtfy = gaa*ay;
      const SIMD_flt_t dtfz = gaa*az;
      const SIMD_flt_t dtgx = fga*ax - hgb*bx;
      const SIMD_flt_t dtgy = fga*ay - hgb*by;
      const SIMD_flt_t dtgz = fga*az - hgb*bz;
      const SIMD_flt_t dthx = gbb*bx;
      const SIMD_flt_t dthy = gbb*by;
      const SIMD_flt_t dthz = gbb*bz;

      const SIMD_flt_t df = -tk * df1;

      const SIMD_flt_t sx2 = df*dtgx;
      const SIMD_flt_t sy2 = df*dtgy;
      const SIMD_flt_t sz2 = df*dtgz;

      SIMD_flt_t f1x = df*dtfx;
      SIMD_flt_t f1y = df*dtfy;
      SIMD_flt_t f1z = df*dtfz;

      SIMD_flt_t f2x = sx2 - f1x;
      SIMD_flt_t f2y = sy2 - f1y;
      SIMD_flt_t f2z = sz2 - f1z;

      SIMD_flt_t f4x = df*dthx;
      SIMD_flt_t f4y = df*dthy;
      SIMD_flt_t f4z = df*dthz;

      SIMD_flt_t f3x = -sx2 - f4x;
      SIMD_flt_t f3y = -sy2 - f4y;
      SIMD_flt_t f3z = -sz2 - f4z;

      SIMD_flt_t qdeng;
      if (EFLAG || VFLAG) {
        SIMD_flt_t ev_pre;
        if (NEWTON_BOND) ev_pre = one;
        else {
          ev_pre = szero;
          const SIMD_flt_t quarter = SIMD_set((flt_t)0.25);
          ev_pre = SIMD_add(ev_pre, i1 < simd_nlocals4, ev_pre, quarter);
          ev_pre = SIMD_add(ev_pre, i2 < simd_nlocals4, ev_pre, quarter);
          ev_pre = SIMD_add(ev_pre, i3 < simd_nlocals4, ev_pre, quarter);
          ev_pre = SIMD_add(ev_pre, i4 < simd_nlocals4, ev_pre, quarter);
        }
        SIMD_zero_masked(nmask, ev_pre);
        if (EFLAG) {
          const SIMD_flt_t deng = tk * p;
          sedihedral = SIMD_ev_add(sedihedral, ev_pre * deng);
          if (eatom) {
            qdeng = deng * SIMD_set((flt_t)0.25);
            SIMD_mask newton_mask;
            if (NEWTON_BOND) newton_mask = nmask;
            if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i2, simd_nlocals4);
            SIMD_flt_t ieng = qdeng;
            SIMD_jeng_update(newton_mask, featom, i2, ieng);
            ieng = qdeng;
            if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i3, simd_nlocals4);
            SIMD_jeng_update(newton_mask, featom, i3, ieng);
          }
        }
        if (VFLAG && vflag) {
          sv0 = SIMD_ev_add(sv0, ev_pre*(vb1x*f1x-vb2xm*f3x+(vb3x-vb2xm)*f4x));
          sv1 = SIMD_ev_add(sv1, ev_pre*(vb1y*f1y-vb2ym*f3y+(vb3y-vb2ym)*f4y));
          sv2 = SIMD_ev_add(sv2, ev_pre*(vb1z*f1z-vb2zm*f3z+(vb3z-vb2zm)*f4z));
          sv3 = SIMD_ev_add(sv3, ev_pre*(vb1x*f1y-vb2xm*f3y+(vb3x-vb2xm)*f4y));
          sv4 = SIMD_ev_add(sv4, ev_pre*(vb1x*f1z-vb2xm*f3z+(vb3x-vb2xm)*f4z));
          sv5 = SIMD_ev_add(sv5, ev_pre*(vb1y*f1z-vb2ym*f3z+(vb3y-vb2ym)*f4z));
        }
      }

      SIMD_mask newton_mask;
      if (NEWTON_BOND) newton_mask = nmask;
      if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i2, simd_nlocals4);
      SIMD_safe_jforce(newton_mask, pforce, i2, f2x, f2y, f2z);
      if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i3, simd_nlocals4);
      SIMD_safe_jforce(newton_mask, pforce, i3, f3x, f3y, f3z);

      // 1-4 LJ and Coulomb interactions
      // tally energy/virial in pair, using newton_bond as newton flag

      const SIMD_flt_t rsq = delx*delx + dely*dely + delz*delz;
      const SIMD_flt_t r2inv = SIMD_rcpz(nmask, rsq);
      const SIMD_flt_t r6inv = r2inv*r2inv*r2inv;

      const SIMD_flt_t simd_qqrd2e = SIMD_set(qqrd2e);
      SIMD_flt_t forcecoul;
      if (implicit) forcecoul = simd_qqrd2e * q1 * q4 * r2inv;
      else forcecoul = simd_qqrd2e * q1 * q4 * SIMD_sqrt(r2inv);

      const SIMD_int ijtype = (itype * ntypes + jtype) << 2;
      const SIMD_flt_t lj1 = SIMD_gather(nmask, plj1, ijtype);
      const SIMD_flt_t lj2 = SIMD_gather(nmask, plj2, ijtype);
      const SIMD_flt_t forcelj = r6inv * (lj1 * r6inv - lj2);
      const SIMD_flt_t fpair = tweight * (forcelj + forcecoul) * r2inv;

      f1x = f1x + delx * fpair;
      f1y = f1y + dely * fpair;
      f1z = f1z + delz * fpair;
      f4x = f4x - delx * fpair;
      f4y = f4y - dely * fpair;
      f4z = f4z - delz * fpair;

      if (EFLAG || VFLAG) {
        SIMD_flt_t ev_pre;
        if (NEWTON_BOND) ev_pre = one;
        else {
          ev_pre = szero;
          const SIMD_flt_t half = SIMD_set((flt_t)0.5);
          ev_pre = SIMD_add(ev_pre, i1 < simd_nlocals4,ev_pre,half);
          ev_pre = SIMD_add(ev_pre, i4 < simd_nlocals4,ev_pre,half);
        }
        SIMD_zero_masked(nmask, ev_pre);

        if (EFLAG) {
          const SIMD_flt_t ecoul = tweight * forcecoul;
          const SIMD_flt_t lj3 = SIMD_gather(nmask, plj3, ijtype);
          const SIMD_flt_t lj4 = SIMD_gather(nmask, plj4, ijtype);
          SIMD_flt_t evdwl = tweight * r6inv * (lj3 * r6inv - lj4);
          secoul = SIMD_ev_add(secoul, ev_pre * ecoul);
          sevdwl = SIMD_ev_add(sevdwl, ev_pre * evdwl);
          if (eatom) {
            const SIMD_flt_t half = SIMD_set((flt_t)0.5);
            evdwl = evdwl * half;
            evdwl = evdwl + half * ecoul + qdeng;

            if (NEWTON_BOND) newton_mask = nmask;
            if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i1, simd_nlocals4);
            SIMD_flt_t ieng = evdwl;
            SIMD_jeng_update(newton_mask, featom, i1, ieng);
            ieng = evdwl;
            if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i4, simd_nlocals4);
            SIMD_jeng_update(newton_mask, featom, i4, ieng);
          }
        }
        if (VFLAG && vflag) {
          spv0 = SIMD_ev_add(spv0, ev_pre * delx * delx * fpair);
          spv1 = SIMD_ev_add(spv1, ev_pre * dely * dely * fpair);
          spv2 = SIMD_ev_add(spv2, ev_pre * delz * delz * fpair);
          spv3 = SIMD_ev_add(spv3, ev_pre * delx * dely * fpair);
          spv4 = SIMD_ev_add(spv4, ev_pre * delx * delz * fpair);
          spv5 = SIMD_ev_add(spv5, ev_pre * dely * delz * fpair);
        }
      }

      if (NEWTON_BOND) newton_mask = nmask;
      if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i1, simd_nlocals4);
      SIMD_safe_jforce(newton_mask, pforce, i1, f1x, f1y, f1z);
      if (!NEWTON_BOND) newton_mask = SIMD_lt(nmask, i4, simd_nlocals4);
      SIMD_safe_jforce(newton_mask, pforce, i4, f4x, f4y, f4z);
    } // for n

    if (EFLAG) {
      oedihedral += SIMD_sum(sedihedral);
      oecoul += SIMD_sum(secoul);
      oevdwl += SIMD_sum(sevdwl);
    }
    if (VFLAG && vflag) {
      ov0 += SIMD_sum(sv0);
      ov1 += SIMD_sum(sv1);
      ov2 += SIMD_sum(sv2);
      ov3 += SIMD_sum(sv3);
      ov4 += SIMD_sum(sv4);
      ov5 += SIMD_sum(sv5);
      opv0 += SIMD_sum(spv0);
      opv1 += SIMD_sum(spv1);
      opv2 += SIMD_sum(spv2);
      opv3 += SIMD_sum(spv3);
      opv4 += SIMD_sum(spv4);
      opv5 += SIMD_sum(spv5);
    }
  } // omp parallel

  if (EFLAG) {
    energy += oedihedral;
    force->pair->eng_vdwl += oevdwl;
    force->pair->eng_coul += oecoul;
  }
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
    force->pair->virial[0] += opv0;
    force->pair->virial[1] += opv1;
    force->pair->virial[2] += opv2;
    force->pair->virial[3] += opv3;
    force->pair->virial[4] += opv4;
    force->pair->virial[5] += opv5;
  }

  fix->set_reduce_flag();
}

#endif

/* ---------------------------------------------------------------------- */

void DihedralCharmmIntel::init_style()
{
  DihedralCharmm::init_style();

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  #ifdef _LMP_INTEL_OFFLOAD
  _use_base = 0;
  if (fix->offload_balance() != 0.0) {
    _use_base = 1;
    return;
  }
  #endif

  fix->bond_init_check();

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void DihedralCharmmIntel::pack_force_const(ForceConst<flt_t> &fc,
                                           IntelBuffers<flt_t,acc_t> *buffers)
{

  const int tp1 = atom->ntypes + 1;
  const int bp1 = atom->ndihedraltypes + 1;
  fc.set_ntypes(tp1,bp1,memory);
  buffers->set_ntypes(tp1);

  if (weightflag) {
    for (int i = 1; i < tp1; i++) {
      for (int j = 1; j < tp1; j++) {
        fc.ljp[i][j].lj1 = lj14_1[i][j];
        fc.ljp[i][j].lj2 = lj14_2[i][j];
        fc.ljp[i][j].lj3 = lj14_3[i][j];
        fc.ljp[i][j].lj4 = lj14_4[i][j];
      }
    }
  }

  for (int i = 1; i < bp1; i++) {
    fc.bp[i].multiplicity = multiplicity[i];
    fc.bp[i].cos_shift = cos_shift[i];
    fc.bp[i].sin_shift = sin_shift[i];
    fc.bp[i].k = k[i];
    fc.weight[i] = weight[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void DihedralCharmmIntel::ForceConst<flt_t>::set_ntypes(const int npairtypes,
                                                        const int nbondtypes,
                                                        Memory *memory) {
  if (npairtypes != _npairtypes) {
    if (_npairtypes > 0)
      _memory->destroy(ljp);
    if (npairtypes > 0)
      memory->create(ljp,npairtypes,npairtypes,"fc.ljp");
  }

  if (nbondtypes != _nbondtypes) {
    if (_nbondtypes > 0) {
      _memory->destroy(bp);
      _memory->destroy(weight);
    }

    if (nbondtypes > 0) {
      _memory->create(bp,nbondtypes,"dihedralcharmmintel.bp");
      _memory->create(weight,nbondtypes,"dihedralcharmmintel.weight");
    }
  }
  _npairtypes = npairtypes;
  _nbondtypes = nbondtypes;
  _memory = memory;
}
