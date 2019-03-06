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

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include "improper_cvff_intel.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "suffix.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
#define SMALL2     (flt_t)0.000001
#define INVSMALL   (flt_t)1000.0
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

ImproperCvffIntel::ImproperCvffIntel(LAMMPS *lmp) :
  ImproperCvff(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

ImproperCvffIntel::~ImproperCvffIntel()
{
}

/* ---------------------------------------------------------------------- */

void ImproperCvffIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    ImproperCvff::compute(eflag, vflag);
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
void ImproperCvffIntel::compute(int eflag, int vflag,
                                    IntelBuffers<flt_t,acc_t> *buffers,
                                    const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

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
void ImproperCvffIntel::eval(const int vflag,
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
  #pragma omp parallel default(none) \
    shared(f_start,f_stride,fc) \
    reduction(+:oeimproper,ov0,ov1,ov2,ov3,ov4,ov5)
  #endif
  {
    int nfrom, npl, nto, tid;
    #ifdef LMP_INTEL_USE_SIMDOFF_FIX
    IP_PRE_omp_range_id(nfrom, nto, tid, inum, nthreads);
    #else
    IP_PRE_omp_stride_id(nfrom, npl, nto, tid, inum, nthreads);
    #endif

    FORCE_T * _noalias const f = f_start + (tid * f_stride);
    if (fix->need_zero(tid))
      memset(f, 0, f_stride * sizeof(FORCE_T));

    const int5_t * _noalias const improperlist =
      (int5_t *) neighbor->improperlist[0];

    #ifdef LMP_INTEL_USE_SIMDOFF_FIX
    acc_t seimproper, sv0, sv1, sv2, sv3, sv4, sv5;
    if (EFLAG) seimproper = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
    }
    #pragma simd reduction(+:seimproper, sv0, sv1, sv2, sv3, sv4, sv5)
    for (int n = nfrom; n < nto; n++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = improperlist[n].a;
      const int i2 = improperlist[n].b;
      const int i3 = improperlist[n].c;
      const int i4 = improperlist[n].d;
      const int type = improperlist[n].t;

      // geometry of 4-body

      const flt_t vb1x = x[i1].x - x[i2].x;
      const flt_t vb1y = x[i1].y - x[i2].y;
      const flt_t vb1z = x[i1].z - x[i2].z;

      const flt_t vb2xm = x[i2].x - x[i3].x;
      const flt_t vb2ym = x[i2].y - x[i3].y;
      const flt_t vb2zm = x[i2].z - x[i3].z;

      const flt_t vb3x = x[i4].x - x[i3].x;
      const flt_t vb3y = x[i4].y - x[i3].y;
      const flt_t vb3z = x[i4].z - x[i3].z;

      // 1st and 2nd angle

      const flt_t b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
      const flt_t rb1 = (flt_t)1.0 / sqrt(b1mag2);
      const flt_t sb1 = (flt_t)1.0 / b1mag2;

      const flt_t b2mag2 = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      const flt_t rb2 = (flt_t)1.0 / sqrt(b2mag2);
      const flt_t sb2 = (flt_t)1.0 / b2mag2;

      const flt_t b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
      const flt_t rb3 = (flt_t)1.0 / sqrt(b3mag2);
      const flt_t sb3 = (flt_t)1.0 / b3mag2;

      const flt_t c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * rb1 * rb3;

      flt_t ctmp = -vb1x*vb2xm - vb1y*vb2ym - vb1z*vb2zm;
      const flt_t r12c1 = rb1 * rb2;
      const flt_t c1mag = ctmp * r12c1;

      ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z;
      const flt_t r12c2 = rb2 * rb3;
      const flt_t c2mag = ctmp * r12c2;

      // cos and sin of 2 angles and final c

      const flt_t sd1 = (flt_t)1.0 - c1mag * c1mag;
      flt_t sc1 = (flt_t)1.0/sqrt(sd1);
      if (sd1 < SMALL2) sc1 = INVSMALL;

      const flt_t sd2 = (flt_t)1.0 - c2mag * c2mag;
      flt_t sc2 = (flt_t)1.0/sqrt(sd2);
      if (sc2 < SMALL2) sc2 = INVSMALL;

      const flt_t s1 = sc1 * sc1;
      const flt_t s2 = sc2 * sc2;
      flt_t s12 = sc1 * sc2;
      flt_t c = (c0 + c1mag*c2mag) * s12;

      // error check
      #ifndef LMP_INTEL_USE_SIMDOFF_FIX
      if (c > PTOLERANCE || c < MTOLERANCE) {
        int me;
        MPI_Comm_rank(world,&me);
        if (screen) {
          char str[128];
          sprintf(str,"Improper problem: %d " BIGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT " "
                  TAGINT_FORMAT " " TAGINT_FORMAT,
                  me,update->ntimestep,
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

      // force & energy
      // p = 1 + cos(n*phi) for d = 1
      // p = 1 - cos(n*phi) for d = -1
      // pd = dp/dc / 2

      const int m = fc.fc[type].multiplicity;

      flt_t p, pd;
      #ifdef LMP_INTEL_USE_SIMDOFF_FIX
      #pragma simdoff
      #endif
      {
        if (m == 2) {
          p = (flt_t)2.0*c*c;
          pd = (flt_t)2.0*c;
        } else if (m == 3) {
          const flt_t rc2 = c*c;
          p = ((flt_t)4.0*rc2-(flt_t)3.0)*c + (flt_t)1.0;
          pd = (flt_t)6.0*rc2 - (flt_t)1.5;
        } else if (m == 4) {
          const flt_t rc2 = c*c;
          p = (flt_t)8.0*(rc2-1)*rc2 + (flt_t)2.0;
          pd = ((flt_t)16.0*rc2-(flt_t)8.0)*c;
        } else if (m == 6) {
          const flt_t rc2 = c*c;
          p = (((flt_t)32.0*rc2-(flt_t)48.0)*rc2 + (flt_t)18.0)*rc2;
          pd = ((flt_t)96.0*(rc2-(flt_t)1.0)*rc2 + (flt_t)18.0)*c;
        } else if (m == 1) {
          p = c + (flt_t)1.0;
          pd = (flt_t)0.5;
        } else if (m == 5) {
          const flt_t rc2 = c*c;
          p = (((flt_t)16.0*rc2-(flt_t)20.0)*rc2 + (flt_t)5.0)*c + (flt_t)1.0;
          pd = ((flt_t)40.0*rc2-(flt_t)30.0)*rc2 + (flt_t)2.5;
        } else if (m == 0) {
          p = (flt_t)2.0;
          pd = (flt_t)0.0;
        }
      }

      if (fc.fc[type].sign == -1) {
        p = (flt_t)2.0 - p;
        pd = -pd;
      }

      flt_t eimproper;
      if (EFLAG) eimproper = fc.fc[type].k * p;

      const flt_t a = (flt_t)2.0 * fc.fc[type].k * pd;
      c = c * a;
      s12 = s12 * a;
      const flt_t a11 = c*sb1*s1;
      const flt_t a22 = -sb2*((flt_t)2.0*c0*s12 - c*(s1+s2));
      const flt_t a33 = c*sb3*s2;
      const flt_t a12 = -r12c1*(c1mag*c*s1 + c2mag*s12);
      const flt_t a13 = -rb1*rb3*s12;
      const flt_t a23 = r12c2*(c2mag*c*s2 + c1mag*s12);

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

      // apply force to each of 4 atoms

      #ifdef LMP_INTEL_USE_SIMDOFF_FIX
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

      if (EFLAG || VFLAG) {
        #ifdef LMP_INTEL_USE_SIMDOFF_FIX
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, eimproper, i1, i2,
                              i3, i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y,
                              f4z, vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm,
                              vb3x, vb3y, vb3z, seimproper, f, NEWTON_BOND,
                              nlocal, sv0, sv1, sv2, sv3, sv4, sv5);
        #else
        IP_PRE_ev_tally_dihed(EFLAG, VFLAG, eatom, vflag, eimproper, i1, i2,
                              i3, i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x, f4y,
                              f4z, vb1x, vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm,
                              vb3x, vb3y, vb3z, oeimproper, f, NEWTON_BOND,
                              nlocal, ov0, ov1, ov2, ov3, ov4, ov5);
        #endif
      }
    } // for n
    #ifdef LMP_INTEL_USE_SIMDOFF_FIX
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

void ImproperCvffIntel::init_style()
{
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
void ImproperCvffIntel::pack_force_const(ForceConst<flt_t> &fc,
                                         IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int bp1 = atom->nimpropertypes + 1;
  fc.set_ntypes(bp1,memory);

  for (int i = 1; i < bp1; i++) {
    fc.fc[i].k = k[i];
    fc.fc[i].sign = sign[i];
    fc.fc[i].multiplicity = multiplicity[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void ImproperCvffIntel::ForceConst<flt_t>::set_ntypes(const int nimproper,
                                                          Memory *memory) {
  if (nimproper != _nimpropertypes) {
    if (_nimpropertypes > 0)
      _memory->destroy(fc);

    if (nimproper > 0)
      _memory->create(fc,nimproper,"improperharmonicintel.fc");
  }
  _nimpropertypes = nimproper;
  _memory = memory;
}
