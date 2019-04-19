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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include "angle_harmonic_intel.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "modify.h"
#include "math_const.h"
#include "memory.h"
#include "suffix.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL2 (flt_t)0.000001
#define INVSMALL (flt_t)1000.0
typedef struct { int a,b,c,t;  } int4_t;

/* ---------------------------------------------------------------------- */

AngleHarmonicIntel::AngleHarmonicIntel(LAMMPS *lmp) : AngleHarmonic(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

AngleHarmonicIntel::~AngleHarmonicIntel()
{
}

/* ---------------------------------------------------------------------- */

void AngleHarmonicIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    AngleHarmonic::compute(eflag, vflag);
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
void AngleHarmonicIntel::compute(int eflag, int vflag,
                               IntelBuffers<flt_t,acc_t> *buffers,
                               const ForceConst<flt_t> &fc)
{
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"USER-INTEL package does not support per-atom stress");

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

/* ---------------------------------------------------------------------- */

template <int EFLAG, int VFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void AngleHarmonicIntel::eval(const int vflag,
                            IntelBuffers<flt_t,acc_t> *buffers,
                            const ForceConst<flt_t> &fc)

{
  const int inum = neighbor->nanglelist;
  if (inum == 0) return;

  ATOM_T * _noalias const x = buffers->get_x(0);
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

  acc_t oeangle, ov0, ov1, ov2, ov3, ov4, ov5;
  if (EFLAG) oeangle = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(f_start,f_stride,fc) \
    reduction(+:oeangle,ov0,ov1,ov2,ov3,ov4,ov5)
  #endif
  {
    int nfrom, npl, nto, tid;
    #ifdef LMP_INTEL_USE_SIMDOFF
    IP_PRE_omp_range_id(nfrom, nto, tid, inum, nthreads);
    #else
    IP_PRE_omp_stride_id(nfrom, npl, nto, tid, inum, nthreads);
    #endif

    FORCE_T * _noalias const f = f_start + (tid * f_stride);
    if (fix->need_zero(tid))
      memset(f, 0, f_stride * sizeof(FORCE_T));

    const int4_t * _noalias const anglelist =
      (int4_t *) neighbor->anglelist[0];

    #ifdef LMP_INTEL_USE_SIMDOFF
    acc_t seangle, sv0, sv1, sv2, sv3, sv4, sv5;
    if (EFLAG) seangle = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
    }
    #pragma simd reduction(+:seangle, sv0, sv1, sv2, sv3, sv4, sv5)
    for (int n = nfrom; n < nto; n ++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = anglelist[n].a;
      const int i2 = anglelist[n].b;
      const int i3 = anglelist[n].c;
      const int type = anglelist[n].t;

      // 1st bond

      const flt_t delx1 = x[i1].x - x[i2].x;
      const flt_t dely1 = x[i1].y - x[i2].y;
      const flt_t delz1 = x[i1].z - x[i2].z;

      const flt_t rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      const flt_t r1 = (flt_t)1.0/sqrt(rsq1);

      // 2nd bond

      const flt_t delx2 = x[i3].x - x[i2].x;
      const flt_t dely2 = x[i3].y - x[i2].y;
      const flt_t delz2 = x[i3].z - x[i2].z;

      const flt_t rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      const flt_t r2 = (flt_t)1.0/sqrt(rsq2);

      // angle (cos and sin)

      flt_t c = delx1*delx2 + dely1*dely2 + delz1*delz2;
      const flt_t r1r2 = r1 * r2;
      c *= r1r2;

      if (c > (flt_t)1.0) c = (flt_t)1.0;
      if (c < (flt_t)-1.0) c = (flt_t)-1.0;

      const flt_t sd = (flt_t)1.0 - c * c;
      flt_t s = (flt_t)1.0/sqrt(sd);
      if (sd < SMALL2) s = INVSMALL;

      // harmonic force & energy

      const flt_t dtheta = acos(c) - fc.fc[type].theta0;
      const flt_t tk = fc.fc[type].k * dtheta;

      flt_t eangle;
      if (EFLAG) eangle = tk*dtheta;

      const flt_t a = (flt_t)-2.0 * tk * s;
      const flt_t ac = a*c;
      const flt_t a11 = ac / rsq1;
      const flt_t a12 = -a * (r1r2);
      const flt_t a22 = ac / rsq2;

      const flt_t f1x = a11*delx1 + a12*delx2;
      const flt_t f1y = a11*dely1 + a12*dely2;
      const flt_t f1z = a11*delz1 + a12*delz2;

      const flt_t f3x = a22*delx2 + a12*delx1;
      const flt_t f3y = a22*dely2 + a12*dely1;
      const flt_t f3z = a22*delz2 + a12*delz1;

      // apply force to each of 3 atoms

      #ifdef LMP_INTEL_USE_SIMDOFF
      #pragma simdoff
      #endif
      {
        if (NEWTON_BOND || i1 < nlocal) {
          f[i1].x += f1x;
          f[i1].y += f1y;
          f[i1].z += f1z;
        }

        if (NEWTON_BOND || i2 < nlocal) {
          f[i2].x -= f1x + f3x;
          f[i2].y -= f1y + f3y;
          f[i2].z -= f1z + f3z;
        }

        if (NEWTON_BOND || i3 < nlocal) {
          f[i3].x += f3x;
          f[i3].y += f3y;
          f[i3].z += f3z;
        }
      }

      if (EFLAG || VFLAG) {
        #ifdef LMP_INTEL_USE_SIMDOFF
        IP_PRE_ev_tally_angle(EFLAG, VFLAG, eatom, vflag, eangle, i1, i2, i3,
                              f1x, f1y, f1z, f3x, f3y, f3z, delx1, dely1,
                              delz1, delx2, dely2, delz2, seangle, f,
                              NEWTON_BOND, nlocal, sv0, sv1, sv2, sv3, sv4,
                              sv5);
        #else
        IP_PRE_ev_tally_angle(EFLAG, VFLAG, eatom, vflag, eangle, i1, i2, i3,
                              f1x, f1y, f1z, f3x, f3y, f3z, delx1, dely1,
                              delz1, delx2, dely2, delz2, oeangle, f,
                              NEWTON_BOND, nlocal, ov0, ov1, ov2, ov3, ov4,
                              ov5);
        #endif
      }
    } // for n
    #ifdef LMP_INTEL_USE_SIMDOFF
    if (EFLAG) oeangle += seangle;
    if (VFLAG && vflag) {
        ov0 += sv0; ov1 += sv1; ov2 += sv2;
        ov3 += sv3; ov4 += sv4; ov5 += sv5;
    }
    #endif
  } // omp parallel

  if (EFLAG) energy += oeangle;
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
  }

  fix->set_reduce_flag();
}

/* ---------------------------------------------------------------------- */

void AngleHarmonicIntel::init_style()
{
  AngleHarmonic::init_style();

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
void AngleHarmonicIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int bp1 = atom->nangletypes + 1;
  fc.set_ntypes(bp1,memory);

  for (int i = 1; i < bp1; i++) {
    fc.fc[i].k = k[i];
    fc.fc[i].theta0 = theta0[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void AngleHarmonicIntel::ForceConst<flt_t>::set_ntypes(const int nangletypes,
                                                     Memory *memory) {
  if (nangletypes != _nangletypes) {
    if (_nangletypes > 0)
      _memory->destroy(fc);

    if (nangletypes > 0)
      _memory->create(fc,nangletypes,"anglecharmmintel.fc");
  }
  _nangletypes = nangletypes;
  _memory = memory;
}
