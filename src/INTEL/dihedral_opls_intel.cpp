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

#include "dihedral_opls_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"

using namespace LAMMPS_NS;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
#define SMALL2     (flt_t)0.000001
#define INVSMALL   (flt_t)1000.0
#define SMALLER2   (flt_t)0.0000000001
#define INVSMALLER (flt_t)100000.0
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

DihedralOPLSIntel::DihedralOPLSIntel(class LAMMPS *lmp)
  : DihedralOPLS(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void DihedralOPLSIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    DihedralOPLS::compute(eflag, vflag);
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
void DihedralOPLSIntel::compute(int eflag, int vflag,
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

template <int EFLAG, int VFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void DihedralOPLSIntel::eval(const int vflag,
                               IntelBuffers<flt_t,acc_t> *buffers,
                               const ForceConst<flt_t> &fc)

{
  const int inum = neighbor->ndihedrallist;
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

  acc_t oedihedral, ov0, ov1, ov2, ov3, ov4, ov5;
  if (EFLAG) oedihedral = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(f_start,f_stride,fc)           \
    reduction(+:oedihedral,ov0,ov1,ov2,ov3,ov4,ov5)
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

    const int5_t * _noalias const dihedrallist =
      (int5_t *) neighbor->dihedrallist[0];

    #ifdef LMP_INTEL_USE_SIMDOFF
    acc_t sedihedral, sv0, sv1, sv2, sv3, sv4, sv5;
    if (EFLAG) sedihedral = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
    }
#if defined(USE_OMP_SIMD)
    #pragma omp simd reduction(+:sedihedral, sv0, sv1, sv2, sv3, sv4, sv5)
#else
    #pragma simd reduction(+:sedihedral, sv0, sv1, sv2, sv3, sv4, sv5)
#endif
    for (int n = nfrom; n < nto; n ++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = IP_PRE_dword_index(dihedrallist[n].a);
      const int i2 = IP_PRE_dword_index(dihedrallist[n].b);
      const int i3 = IP_PRE_dword_index(dihedrallist[n].c);
      const int i4 = IP_PRE_dword_index(dihedrallist[n].d);
      const int type = IP_PRE_dword_index(dihedrallist[n].t);

      // 1st bond

      const flt_t vb1x = x[i1].x - x[i2].x;
      const flt_t vb1y = x[i1].y - x[i2].y;
      const flt_t vb1z = x[i1].z - x[i2].z;

      // 2nd bond

      const flt_t vb2xm = x[i2].x - x[i3].x;
      const flt_t vb2ym = x[i2].y - x[i3].y;
      const flt_t vb2zm = x[i2].z - x[i3].z;

      // 3rd bond

      const flt_t vb3x = x[i4].x - x[i3].x;
      const flt_t vb3y = x[i4].y - x[i3].y;
      const flt_t vb3z = x[i4].z - x[i3].z;

      // c0 calculation
      // 1st and 2nd angle

      const flt_t b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
      const flt_t rb1 = (flt_t)1.0 / std::sqrt(b1mag2);
      const flt_t sb1 = (flt_t)1.0 / b1mag2;

      const flt_t b2mag2 = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      const flt_t rb2 = (flt_t)1.0 / std::sqrt(b2mag2);
      const flt_t sb2 = (flt_t)1.0 / b2mag2;

      const flt_t b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
      const flt_t rb3 = (flt_t)1.0 / std::sqrt(b3mag2);
      const flt_t sb3 = (flt_t)1.0 / b3mag2;

      const flt_t c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3;

      flt_t ctmp = -vb1x*vb2xm - vb1y*vb2ym - vb1z*vb2zm;
      const flt_t r12c1 =  rb1 * rb2;
      const flt_t c1mag = ctmp * r12c1;

      ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
      const flt_t r12c2 =  rb2 * rb3;
      const flt_t c2mag = ctmp * r12c2;

      // cos and sin of 2 angles and final c

      flt_t sin2 = MAX((flt_t)1.0 - c1mag*c1mag,(flt_t)0.0);
      flt_t sc1 = (flt_t)1.0/std::sqrt(sin2);
      if (sin2 < SMALL2) sc1 = INVSMALL;

      sin2 = MAX((flt_t)1.0 - c2mag*c2mag,(flt_t)0.0);
      flt_t sc2 = (flt_t)1.0/std::sqrt(sin2);
      if (sin2 < SMALL2) sc2 = INVSMALL;

      const flt_t s1 = sc1 * sc1;
      const flt_t s2 = sc2 * sc2;
      flt_t s12 = sc1 * sc2;
      flt_t c = (c0 + c1mag*c2mag) * s12;

      const flt_t cx = vb1z*vb2ym - vb1y*vb2zm;
      const flt_t cy = vb1x*vb2zm - vb1z*vb2xm;
      const flt_t cz = vb1y*vb2xm - vb1x*vb2ym;
      const flt_t cmag = (flt_t)1.0/std::sqrt(cx*cx + cy*cy + cz*cz);
      const flt_t dx = (cx*vb3x + cy*vb3y + cz*vb3z)*cmag*rb3;

      // error check
      #ifndef LMP_INTEL_USE_SIMDOFF
      if (c > PTOLERANCE || c < MTOLERANCE)
        problem(FLERR, i1, i2, i3, i4);
      #endif

      if (c > (flt_t)1.0) c = (flt_t)1.0;
      if (c < (flt_t)-1.0) c = (flt_t)-1.0;

      // force & energy
      // p = sum (i=1,4) k_i * (1 + (-1)**(i+1)*cos(i*phi) )
      // pd = dp/dc

      const flt_t cossq = c * c;
      const flt_t sinsq = (flt_t)1.0 - cossq;
      flt_t siinv = (flt_t)1.0/std::sqrt(sinsq);
      if (sinsq < SMALLER2 ) siinv = INVSMALLER;
      if (dx < (flt_t)0.0) siinv = -siinv;

      const flt_t cos_2phi = cossq - sinsq;
      const flt_t sin_2phim = (flt_t)2.0 * c;
      const flt_t cos_3phi = (flt_t)2.0 * c * cos_2phi - c;
      const flt_t sin_3phim = (flt_t)2.0 * cos_2phi + (flt_t)1.0;
      const flt_t cos_4phi = (flt_t)2.0 * cos_2phi * cos_2phi - (flt_t)1.0;
      const flt_t sin_4phim = (flt_t)2.0 * cos_2phi * sin_2phim;

      flt_t p, pd;
      p = fc.fc[type].k1*((flt_t)1.0 + c) +
          fc.fc[type].k2*((flt_t)1.0 - cos_2phi) +
          fc.fc[type].k3*((flt_t)1.0 + cos_3phi) +
          fc.fc[type].k4*((flt_t)1.0 - cos_4phi) ;
      pd = fc.fc[type].k1 -
           (flt_t)2.0 * fc.fc[type].k2 * sin_2phim +
           (flt_t)3.0 * fc.fc[type].k3 * sin_3phim -
           (flt_t)4.0 * fc.fc[type].k4 * sin_4phim;

      flt_t edihed;
      if (EFLAG) edihed = p;

      const flt_t a = pd;
      c = c * a;
      s12 = s12 * a;
      const flt_t a11 = c*sb1*s1;
      const flt_t a22 = -sb2 * ((flt_t)2.0*c0*s12 - c*(s1+s2));
      const flt_t a33 = c*sb3*s2;
      const flt_t a12 = -r12c1 * (c1mag*c*s1 + c2mag*s12);
      const flt_t a13 = -rb1*rb3*s12;
      const flt_t a23 = r12c2 * (c2mag*c*s2 + c1mag*s12);

      const flt_t sx2  = a12*vb1x - a22*vb2xm + a23*vb3x;
      const flt_t sy2  = a12*vb1y - a22*vb2ym + a23*vb3y;
      const flt_t sz2  = a12*vb1z - a22*vb2zm + a23*vb3z;

      const flt_t f1x = a11*vb1x - a12*vb2xm + a13*vb3x;
      const flt_t f1y = a11*vb1y - a12*vb2ym + a13*vb3y;
      const flt_t f1z = a11*vb1z - a12*vb2zm + a13*vb3z;

      const flt_t f2x = -sx2 - f1x;
      const flt_t f2y = -sy2 - f1y;
      const flt_t f2z = -sz2 - f1z;

      const flt_t f4x = a13*vb1x - a23*vb2xm + a33*vb3x;
      const flt_t f4y = a13*vb1y - a23*vb2ym + a33*vb3y;
      const flt_t f4z = a13*vb1z - a23*vb2zm + a33*vb3z;

      const flt_t f3x = sx2 - f4x;
      const flt_t f3y = sy2 - f4y;
      const flt_t f3z = sz2 - f4z;

      if (EFLAG || VFLAG) {
        #ifdef LMP_INTEL_USE_SIMDOFF
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, edihed, i1, i2, i3,
                              i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y, f4z,
                              vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm, vb3x,
                              vb3y, vb3z, sedihedral, f, NEWTON_BOND, nlocal,
                              sv0, sv1, sv2, sv3, sv4, sv5);
        #else
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, edihed, i1, i2, i3,
                              i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y, f4z,
                              vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm, vb3x,
                              vb3y, vb3z, oedihedral, f, NEWTON_BOND, nlocal,
                              ov0, ov1, ov2, ov3, ov4, ov5);
        #endif
      }

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
    } // for n
    #ifdef LMP_INTEL_USE_SIMDOFF
    if (EFLAG) oedihedral += sedihedral;
    if (VFLAG && vflag) {
        ov0 += sv0; ov1 += sv1; ov2 += sv2;
        ov3 += sv3; ov4 += sv4; ov5 += sv5;
    }
    #endif
  } // omp parallel

  if (EFLAG) energy += oedihedral;
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
  }

  fix->set_reduce_flag();
}

/* ---------------------------------------------------------------------- */

void DihedralOPLSIntel::init_style()
{
  DihedralOPLS::init_style();

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
void DihedralOPLSIntel::pack_force_const(ForceConst<flt_t> &fc,
                                         IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int dp1 = atom->ndihedraltypes + 1;
  fc.set_ntypes(dp1,memory);

  for (int i = 1; i < dp1; i++) {
    fc.fc[i].k1 = k1[i];
    fc.fc[i].k2 = k2[i];
    fc.fc[i].k3 = k3[i];
    fc.fc[i].k4 = k4[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void DihedralOPLSIntel::ForceConst<flt_t>::set_ntypes(const int ndihderaltypes,
                                                          Memory *memory) {
  if (memory != nullptr) _memory = memory;
  if (ndihderaltypes != _ndihderaltypes) {
    _memory->destroy(fc);

    if (ndihderaltypes > 0)
      _memory->create(fc,ndihderaltypes,"dihedralcharmmintel.fc");
  }
  _ndihderaltypes = ndihderaltypes;
}
