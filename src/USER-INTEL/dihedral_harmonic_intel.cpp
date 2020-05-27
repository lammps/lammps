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

#include "omp_compat.h"
#include <mpi.h>
#include <cmath>
#include "dihedral_harmonic_intel.h"
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

#include "suffix.h"
using namespace LAMMPS_NS;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

DihedralHarmonicIntel::DihedralHarmonicIntel(class LAMMPS *lmp)
  : DihedralHarmonic(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void DihedralHarmonicIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    DihedralHarmonic::compute(eflag, vflag);
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
void DihedralHarmonicIntel::compute(int eflag, int vflag,
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

template <int EFLAG, int VFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void DihedralHarmonicIntel::eval(const int vflag,
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
    #pragma simd reduction(+:sedihedral, sv0, sv1, sv2, sv3, sv4, sv5)
    for (int n = nfrom; n < nto; n ++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = dihedrallist[n].a;
      const int i2 = dihedrallist[n].b;
      const int i3 = dihedrallist[n].c;
      const int i4 = dihedrallist[n].d;
      const int type = dihedrallist[n].t;

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
      #ifndef LMP_INTEL_USE_SIMDOFF
      if (c > PTOLERANCE || c < MTOLERANCE) {
        int me = comm->me;

        if (screen) {
          char str[128];
          sprintf(str,"Dihedral problem: %d/%d " BIGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT,
                  me,tid,update->ntimestep,
                  atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
          error->warning(FLERR,str,0);
          fprintf(screen,"  1st atom: %d %g %g %g\n",
                  me,x[i1].x,x[i1].y,x[i1].z);
          fprintf(screen,"  2nd atom: %d %g %g %g\n",
                  me,x[i2].x,x[i2].y,x[i2].z);
          fprintf(screen,"  3rd atom: %d %g %g %g\n",
                  me,x[i3].x,x[i3].y,x[i3].z);
          fprintf(screen,"  4th atom: %d %g %g %g\n",
                  me,x[i4].x,x[i4].y,x[i4].z);
        }
      }
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
        #ifdef LMP_INTEL_USE_SIMDOFF
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, deng, i1, i2, i3, i4,
                              f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y, f4z,
                              vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm, vb3x,
                              vb3y, vb3z, sedihedral, f, NEWTON_BOND, nlocal,
                              sv0, sv1, sv2, sv3, sv4, sv5);
        #else
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, deng, i1, i2, i3, i4,
                              f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y, f4z,
                              vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm, vb3x,
                              vb3y, vb3z, oedihedral, f, NEWTON_BOND, nlocal,
                              ov0, ov1, ov2, ov3, ov4, ov5);
        #endif
      }

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

void DihedralHarmonicIntel::init_style()
{
  DihedralHarmonic::init_style();

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
void DihedralHarmonicIntel::pack_force_const(ForceConst<flt_t> &fc,
                                             IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int bp1 = atom->ndihedraltypes + 1;
  fc.set_ntypes(bp1,memory);

  for (int i = 1; i < bp1; i++) {
    fc.bp[i].multiplicity = multiplicity[i];
    fc.bp[i].cos_shift = cos_shift[i];
    fc.bp[i].sin_shift = sin_shift[i];
    fc.bp[i].k = k[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void DihedralHarmonicIntel::ForceConst<flt_t>::set_ntypes(const int nbondtypes,
                                                          Memory *memory) {
  if (nbondtypes != _nbondtypes) {
    if (_nbondtypes > 0)
      _memory->destroy(bp);

    if (nbondtypes > 0)
      _memory->create(bp,nbondtypes,"dihedralcharmmintel.bp");
  }
  _nbondtypes = nbondtypes;
  _memory = memory;
}
