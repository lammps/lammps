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
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */


#include "bond_fene_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using MathConst::MY_CUBEROOT2;

typedef struct { int a,b,t;  } int3_t;

/* ---------------------------------------------------------------------- */

BondFENEIntel::BondFENEIntel(LAMMPS *lmp) : BondFENE(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

BondFENEIntel::~BondFENEIntel()
{
}

/* ---------------------------------------------------------------------- */

void BondFENEIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    BondFENE::compute(eflag, vflag);
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
void BondFENEIntel::compute(int eflag, int vflag,
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
void BondFENEIntel::eval(const int vflag,
                         IntelBuffers<flt_t,acc_t> *buffers,
                         const ForceConst<flt_t> &fc)
{
  const int inum = neighbor->nbondlist;
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

  acc_t oebond, ov0, ov1, ov2, ov3, ov4, ov5;
  if (EFLAG) oebond = (acc_t)0.0;
  if (VFLAG && vflag) {
    ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(f_start,f_stride,fc)           \
    reduction(+:oebond,ov0,ov1,ov2,ov3,ov4,ov5)
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

    const int3_t * _noalias const bondlist =
      (int3_t *) neighbor->bondlist[0];

    #ifdef LMP_INTEL_USE_SIMDOFF
    acc_t sebond, sv0, sv1, sv2, sv3, sv4, sv5;
    if (EFLAG) sebond = (acc_t)0.0;
    if (VFLAG && vflag) {
      sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
    }
#if defined(USE_OMP_SIMD)
    #pragma omp simd reduction(+:sebond, sv0, sv1, sv2, sv3, sv4, sv5)
#else
    #pragma simd reduction(+:sebond, sv0, sv1, sv2, sv3, sv4, sv5)
#endif
    for (int n = nfrom; n < nto; n ++) {
    #else
    for (int n = nfrom; n < nto; n += npl) {
    #endif
      const int i1 = bondlist[n].a;
      const int i2 = bondlist[n].b;
      const int type = bondlist[n].t;

      const flt_t ir0sq = fc.fc[type].ir0sq;
      const flt_t k = fc.fc[type].k;
      const flt_t sigma = fc.fc[type].sigma;
      const flt_t sigmasq = sigma*sigma;
      const flt_t epsilon = fc.fc[type].epsilon;

      const flt_t delx = x[i1].x - x[i2].x;
      const flt_t dely = x[i1].y - x[i2].y;
      const flt_t delz = x[i1].z - x[i2].z;

      const flt_t rsq = delx*delx + dely*dely + delz*delz;
      flt_t rlogarg = (flt_t)1.0 - rsq * ir0sq;
      flt_t irsq = (flt_t)1.0 / rsq;

      // if r -> r0, then rlogarg < 0.0 which is an error
      // issue a warning and reset rlogarg = epsilon
      // if r > 2*r0 something serious is wrong, abort

      if (rlogarg < (flt_t)0.1) {
        error->warning(FLERR,"FENE bond too long: {} {} {} {:.8}",
                       update->ntimestep,atom->tag[i1],atom->tag[i2],sqrt(rsq));
        if (rlogarg <= (flt_t)-3.0) error->one(FLERR,"Bad FENE bond");
        rlogarg = (flt_t)0.1;
      }

      flt_t fbond = -k/rlogarg;

      // force from LJ term

      flt_t sr2,sr6;
      if (rsq < (flt_t)MY_CUBEROOT2*sigmasq) {
        sr2 = sigmasq * irsq;
        sr6 = sr2 * sr2 * sr2;
        fbond += (flt_t)48.0 * epsilon * sr6 * (sr6 - (flt_t)0.5) * irsq;
      }

      // energy

      flt_t ebond;
      if (EFLAG) {
        ebond = (flt_t)-0.5 * k / ir0sq * log(rlogarg);
        if (rsq < (flt_t)MY_CUBEROOT2 * sigmasq)
          ebond += (flt_t)4.0 * epsilon * sr6 * (sr6 - (flt_t)1.0) + epsilon;
      }

      // apply force to each of 2 atoms

      #ifdef LMP_INTEL_USE_SIMDOFF
#if defined(USE_OMP_SIMD)
      #pragma omp ordered simd
#else
      #pragma simdoff
#endif
      #endif
      {
        if (NEWTON_BOND || i1 < nlocal) {
          f[i1].x += delx*fbond;
          f[i1].y += dely*fbond;
          f[i1].z += delz*fbond;
        }

        if (NEWTON_BOND || i2 < nlocal) {
          f[i2].x -= delx*fbond;
          f[i2].y -= dely*fbond;
          f[i2].z -= delz*fbond;
        }
      }

      if (EFLAG || VFLAG) {
        #ifdef LMP_INTEL_USE_SIMDOFF
        IP_PRE_ev_tally_bond(EFLAG, VFLAG, eatom, vflag, ebond, i1, i2, fbond,
                             delx, dely, delz, sebond, f, NEWTON_BOND,
                             nlocal, sv0, sv1, sv2, sv3, sv4, sv5);
        #else
        IP_PRE_ev_tally_bond(EFLAG, VFLAG, eatom, vflag, ebond, i1, i2, fbond,
                             delx, dely, delz, oebond, f, NEWTON_BOND,
                             nlocal, ov0, ov1, ov2, ov3, ov4, ov5);
        #endif
      }
    } // for n
    #ifdef LMP_INTEL_USE_SIMDOFF
    if (EFLAG) oebond += sebond;
    if (VFLAG && vflag) {
       ov0 += sv0; ov1 += sv1; ov2 += sv2;
       ov3 += sv3; ov4 += sv4; ov5 += sv5;
    }
    #endif
  } // omp parallel

  if (EFLAG) energy += oebond;
  if (VFLAG && vflag) {
    virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
    virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
  }

  fix->set_reduce_flag();
}

/* ---------------------------------------------------------------------- */

void BondFENEIntel::init_style()
{
  BondFENE::init_style();

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
void BondFENEIntel::pack_force_const(ForceConst<flt_t> &fc,
                                     IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  const int bp1 = atom->nbondtypes + 1;
  fc.set_ntypes(bp1,memory);

  for (int i = 1; i < bp1; i++) {
    fc.fc[i].k = k[i];
    fc.fc[i].ir0sq = 1.0 / (r0[i] * r0[i]);
    fc.fc[i].sigma = sigma[i];
    fc.fc[i].epsilon = epsilon[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void BondFENEIntel::ForceConst<flt_t>::set_ntypes(const int nbondtypes,
                                                      Memory *memory) {
  if (nbondtypes != _nbondtypes) {
    if (_nbondtypes > 0)
      _memory->destroy(fc);

    if (nbondtypes > 0)
      _memory->create(fc,nbondtypes,"bondfeneintel.fc");
  }
  _nbondtypes = nbondtypes;
  _memory = memory;
}
