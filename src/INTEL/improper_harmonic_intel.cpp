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


#include "improper_harmonic_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
#define SMALL     (flt_t)0.001
#define SMALL2     (flt_t)0.000001
#define INVSMALL   (flt_t)1000.0
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

ImproperHarmonicIntel::ImproperHarmonicIntel(LAMMPS *lmp) :
  ImproperHarmonic(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonicIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    ImproperHarmonic::compute(eflag, vflag);
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
void ImproperHarmonicIntel::compute(int eflag, int vflag,
                                    IntelBuffers<flt_t,acc_t> *buffers,
                                    const ForceConst<flt_t> &fc)
{
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"INTEL package does not support per-atom stress");

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
void ImproperHarmonicIntel::eval(const int vflag,
                                 IntelBuffers<flt_t,acc_t> *buffers,
                                 const ForceConst<flt_t> &fc)
{
  const int inum = neighbor->nimproperlist;
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

  acc_t oeimproper, ov0, ov1, ov2, ov3, ov4, ov5;
  if (EFLAG) oeimproper = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(f_start,f_stride,fc) \
    reduction(+:oeimproper,ov0,ov1,ov2,ov3,ov4,ov5)
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

    const int5_t * _noalias const improperlist =
      (int5_t *) neighbor->improperlist[0];

    #ifdef LMP_INTEL_USE_SIMDOFF
    acc_t seimproper, sv0, sv1, sv2, sv3, sv4, sv5;
    if (EFLAG) seimproper = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
    }
#if defined(USE_OMP_SIMD)
    #pragma omp simd reduction(+:seimproper, sv0, sv1, sv2, sv3, sv4, sv5)
#else
    #pragma simd reduction(+:seimproper, sv0, sv1, sv2, sv3, sv4, sv5)
#endif
    for (int n = nfrom; n < nto; n++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = IP_PRE_dword_index(improperlist[n].a);
      const int i2 = IP_PRE_dword_index(improperlist[n].b);
      const int i3 = IP_PRE_dword_index(improperlist[n].c);
      const int i4 = IP_PRE_dword_index(improperlist[n].d);
      const int type = IP_PRE_dword_index(improperlist[n].t);

      // geometry of 4-body

      const flt_t vb1x = x[i1].x - x[i2].x;
      const flt_t vb1y = x[i1].y - x[i2].y;
      const flt_t vb1z = x[i1].z - x[i2].z;

      const flt_t vb2x = x[i3].x - x[i2].x;
      const flt_t vb2y = x[i3].y - x[i2].y;
      const flt_t vb2z = x[i3].z - x[i2].z;

      const flt_t vb3x = x[i4].x - x[i3].x;
      const flt_t vb3y = x[i4].y - x[i3].y;
      const flt_t vb3z = x[i4].z - x[i3].z;

      flt_t ss1 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
      flt_t ss2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
      flt_t ss3 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;

      const flt_t r1 = (flt_t)1.0 / sqrt(ss1);
      const flt_t r2 = (flt_t)1.0 / sqrt(ss2);
      const flt_t r3 = (flt_t)1.0 / sqrt(ss3);

      ss1 = (flt_t)1.0 / ss1;
      ss2 = (flt_t)1.0 / ss2;
      ss3 = (flt_t)1.0 / ss3;

      // sin and cos of angle

      const flt_t c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
      const flt_t c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
      const flt_t c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

      flt_t s1 = 1.0 - c1*c1;
      if (s1 < SMALL) s1 = SMALL;

      flt_t s2 = (flt_t)1.0 - c2*c2;
      if (s2 < SMALL) s2 = SMALL;

      flt_t s12 = (flt_t)1.0 / sqrt(s1*s2);
      s1 = (flt_t)1.0 / s1;
      s2 = (flt_t)1.0 / s2;
      flt_t c = (c1*c2 + c0) * s12;

      // error check
      #ifndef LMP_INTEL_USE_SIMDOFF
      if (c > PTOLERANCE || c < MTOLERANCE)
        problem(FLERR, i1, i2, i3, i4);
      #endif

      if (c > (flt_t)1.0) c = (flt_t)1.0;
      if (c < (flt_t)-1.0) c = (flt_t)-1.0;

      const flt_t sd = (flt_t)1.0 - c * c;
      flt_t s = (flt_t)1.0 / sqrt(sd);
      if (sd < SMALL2) s = INVSMALL;

      // force & energy

      const flt_t domega = acos(c) - fc.fc[type].chi;
      flt_t a;
      a = fc.fc[type].k * domega;

      flt_t eimproper;
      if (EFLAG) eimproper = a*domega;

      a = -a * (flt_t)2.0 * s;
      c = c * a;
      s12 = s12 * a;
      const flt_t a11 = c*ss1*s1;
      const flt_t a22 = -ss2 * ((flt_t)2.0*c0*s12 - c*(s1+s2));
      const flt_t a33 = c*ss3*s2;
      const flt_t a12 = -r1*r2*(c1*c*s1 + c2*s12);
      const flt_t a13 = -r1*r3*s12;
      const flt_t a23 = r2*r3*(c2*c*s2 + c1*s12);

      const flt_t sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
      const flt_t sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
      const flt_t sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

      const flt_t f1x = a12*vb2x + a13*vb3x + a11*vb1x;
      const flt_t f1y = a12*vb2y + a13*vb3y + a11*vb1y;
      const flt_t f1z = a12*vb2z + a13*vb3z + a11*vb1z;

      const flt_t f2x = -sx2 - f1x;
      const flt_t f2y = -sy2 - f1y;
      const flt_t f2z = -sz2 - f1z;

      const flt_t f4x = a23*vb2x + a33*vb3x + a13*vb1x;
      const flt_t f4y = a23*vb2y + a33*vb3y + a13*vb1y;
      const flt_t f4z = a23*vb2z + a33*vb3z + a13*vb1z;

      const flt_t f3x = sx2 - f4x;
      const flt_t f3y = sy2 - f4y;
      const flt_t f3z = sz2 - f4z;

      // apply force to each of 4 atoms

      #ifdef LMP_INTEL_USE_SIMDOFF
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

        if (NEWTON_BOND || i4 < nlocal) {
          f[i4].x += f4x;
          f[i4].y += f4y;
          f[i4].z += f4z;
        }
      }

      if (EFLAG || VFLAG) {
        #ifdef LMP_INTEL_USE_SIMDOFF
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, eimproper, i1, i2,
                              i3, i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x,
                              f4y, f4z, vb1x, vb1y, vb1z, vb2x, vb2y, vb2z,
                              vb3x, vb3y, vb3z, seimproper, f, NEWTON_BOND,
                              nlocal, sv0, sv1, sv2, sv3, sv4, sv5);
        #else
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, eimproper, i1, i2,
                              i3, i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x,
                              f4y, f4z, vb1x, vb1y, vb1z, vb2x, vb2y, vb2z,
                              vb3x, vb3y, vb3z, oeimproper, f, NEWTON_BOND,
                              nlocal, ov0, ov1, ov2, ov3, ov4, ov5);
        #endif
      }
    } // for n
    #ifdef LMP_INTEL_USE_SIMDOFF
    if (EFLAG) oeimproper += seimproper;
    if (VFLAG && vflag) {
      ov0 += sv0; ov1 += sv1; ov2 += sv2;
      ov3 += sv3; ov4 += sv4; ov5 += sv5;
    }
    #endif
  } // omp parallel
  if (EFLAG) energy += oeimproper;
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
  }

  fix->set_reduce_flag();
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonicIntel::init_style()
{
  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

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
void ImproperHarmonicIntel::pack_force_const(ForceConst<flt_t> &fc,
                                             IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int ip1 = atom->nimpropertypes + 1;
  fc.set_ntypes(ip1,memory);

  for (int i = 1; i < ip1; i++) {
    fc.fc[i].k = k[i];
    fc.fc[i].chi = chi[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void ImproperHarmonicIntel::ForceConst<flt_t>::set_ntypes(const int nimpropertypes,
                                                          Memory *memory) {
  if (memory != nullptr) _memory = memory;
  if (nimpropertypes != _nimpropertypes) {
    _memory->destroy(fc);

    if (nimpropertypes > 0)
      _memory->create(fc,nimpropertypes,"improperharmonicintel.fc");
  }
  _nimpropertypes = nimpropertypes;
}
