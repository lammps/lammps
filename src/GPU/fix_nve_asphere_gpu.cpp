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

#include "fix_nve_asphere_gpu.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include <cmath>
#if (LAL_USE_OMP == 1)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

#define ME_qnormalize(q)                                                \
{                                                                       \
  double norm = 1.0 /                                                   \
    sqrt(q##_w*q##_w + q##_i*q##_i + q##_j*q##_j + q##_k*q##_k);        \
  q##_w *= norm;                                                        \
  q##_i *= norm;                                                        \
  q##_j *= norm;                                                        \
  q##_k *= norm;                                                        \
}

#define ME_mq_to_omega(m, quat, moments_0, moments_1, moments_2, w)     \
{                                                                       \
  double wbody_0, wbody_1, wbody_2;                                     \
  double rot_0, rot_1, rot_2, rot_3, rot_4, rot_5, rot_6, rot_7, rot_8; \
                                                                        \
  double w2 = quat##_w * quat##_w;                                      \
  double i2 = quat##_i * quat##_i;                                      \
  double j2 = quat##_j * quat##_j;                                      \
  double k2 = quat##_k * quat##_k;                                      \
  double twoij = 2.0 * quat##_i * quat##_j;                             \
  double twoik = 2.0 * quat##_i * quat##_k;                             \
  double twojk = 2.0 * quat##_j * quat##_k;                             \
  double twoiw = 2.0 * quat##_i * quat##_w;                             \
  double twojw = 2.0 * quat##_j * quat##_w;                             \
  double twokw = 2.0 * quat##_k * quat##_w;                             \
                                                                        \
  rot##_0 = w2 + i2 - j2 - k2;                                          \
  rot##_1 = twoij - twokw;                                              \
  rot##_2 = twojw + twoik;                                              \
                                                                        \
  rot##_3 = twoij + twokw;                                              \
  rot##_4 = w2 - i2 + j2 - k2;                                          \
  rot##_5 = twojk - twoiw;                                              \
                                                                        \
  rot##_6 = twoik - twojw;                                              \
  rot##_7 = twojk + twoiw;                                              \
  rot##_8 = w2 - i2 - j2 + k2;                                          \
                                                                        \
  wbody_0 = rot##_0*m##_0 + rot##_3*m##_1 + rot##_6*m##_2;              \
  wbody_1 = rot##_1*m##_0 + rot##_4*m##_1 + rot##_7*m##_2;              \
  wbody_2 = rot##_2*m##_0 + rot##_5*m##_1 + rot##_8*m##_2;              \
                                                                        \
  wbody_0 *= moments_0;                                                 \
  wbody_1 *= moments_1;                                                 \
  wbody_2 *= moments_2;                                                 \
                                                                        \
  w##_0 = rot##_0*wbody_0 + rot##_1*wbody_1 + rot##_2*wbody_2;          \
  w##_1 = rot##_3*wbody_0 + rot##_4*wbody_1 + rot##_5*wbody_2;          \
  w##_2 = rot##_6*wbody_0 + rot##_7*wbody_1 + rot##_8*wbody_2;          \
}

#define ME_omega_richardson(dtf,dtq,angmomin,quatin,torque,i0,i1,i2)    \
{                                                                       \
  angmomin[0] += dtf * torque[0];                                       \
  double angmom_0 = angmomin[0];                                        \
  angmomin[1] += dtf * torque[1];                                       \
  double angmom_1 = angmomin[1];                                        \
  angmomin[2] += dtf * torque[2];                                       \
  double angmom_2 = angmomin[2];                                        \
                                                                        \
  double quat_w = quatin[0];                                            \
  double quat_i = quatin[1];                                            \
  double quat_j = quatin[2];                                            \
  double quat_k = quatin[3];                                            \
                                                                        \
  double omega_0, omega_1, omega_2;                                     \
  ME_mq_to_omega(angmom,quat,i0,i1,i2,omega);                           \
                                                                        \
  double wq_0, wq_1, wq_2, wq_3;                                        \
  wq_0 = -omega_0*quat_i - omega_1*quat_j - omega_2*quat_k;             \
  wq_1 = quat_w*omega_0 + omega_1*quat_k - omega_2*quat_j;              \
  wq_2 = quat_w*omega_1 + omega_2*quat_i - omega_0*quat_k;              \
  wq_3 = quat_w*omega_2 + omega_0*quat_j - omega_1*quat_i;              \
                                                                        \
  double qfull_w, qfull_i, qfull_j, qfull_k;                            \
  qfull_w = quat_w + dtq * wq_0;                                        \
  qfull_i = quat_i + dtq * wq_1;                                        \
  qfull_j = quat_j + dtq * wq_2;                                        \
  qfull_k = quat_k + dtq * wq_3;                                        \
  ME_qnormalize(qfull);                                                 \
                                                                        \
  double qhalf_w, qhalf_i, qhalf_j, qhalf_k;                            \
  qhalf_w = quat_w + 0.5*dtq * wq_0;                                    \
  qhalf_i = quat_i + 0.5*dtq * wq_1;                                    \
  qhalf_j = quat_j + 0.5*dtq * wq_2;                                    \
  qhalf_k = quat_k + 0.5*dtq * wq_3;                                    \
  ME_qnormalize(qhalf);                                                 \
                                                                        \
  ME_mq_to_omega(angmom,qhalf,i0,i1,i2,omega);                          \
  wq_0 = -omega_0*qhalf_i - omega_1*qhalf_j - omega_2*qhalf_k;          \
  wq_1 = qhalf_w*omega_0 + omega_1*qhalf_k - omega_2*qhalf_j;           \
  wq_2 = qhalf_w*omega_1 + omega_2*qhalf_i - omega_0*qhalf_k;           \
  wq_3 = qhalf_w*omega_2 + omega_0*qhalf_j - omega_1*qhalf_i;           \
                                                                        \
  qhalf_w += 0.5*dtq * wq_0;                                            \
  qhalf_i += 0.5*dtq * wq_1;                                            \
  qhalf_j += 0.5*dtq * wq_2;                                            \
  qhalf_k += 0.5*dtq * wq_3;                                            \
  ME_qnormalize(qhalf);                                                 \
                                                                        \
  quat_w = 2.0*qhalf_w - qfull_w;                                       \
  quat_i = 2.0*qhalf_i - qfull_i;                                       \
  quat_j = 2.0*qhalf_j - qfull_j;                                       \
  quat_k = 2.0*qhalf_k - qfull_k;                                       \
  ME_qnormalize(quat);                                                  \
                                                                        \
  quatin[0] = quat_w;                                                   \
  quatin[1] = quat_i;                                                   \
  quatin[2] = quat_j;                                                   \
  quatin[3] = quat_k;                                                   \
}

/* ---------------------------------------------------------------------- */

FixNVEAsphereGPU::FixNVEAsphereGPU(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  _dtfm = nullptr;
  _nlocal_max = 0;
  _inertia0 = nullptr;
  _inertia1 = nullptr;
  _inertia2 = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereGPU::init()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec)
    error->all(FLERR,"Compute nve/asphere requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/asphere requires extended particles");

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereGPU::setup(int vflag)
{
  FixNVE::setup(vflag);
  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereGPU::initial_integrate(int /*vflag*/)
{
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double * _noalias const x = atom->x[0];
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  int *mask = atom->mask;

  double **angmom = atom->angmom;
  double **torque = atom->torque;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = comm->nthreads;
    const int tid = omp_get_thread_num();
    const int idelta = nlocal / nthreads + 1;
    const int ifrom = tid * idelta;
    const int ito = MIN(ifrom + idelta, nlocal);
    const int ifrom3 = ifrom * 3;
    const int ito3 = ito * 3;
    #else
    const int tid = 0;
    const int ifrom = 0;
    const int ifrom3 = 0;
    const int ito = nlocal;
    const int ito3 = nlocal * 3;
    #endif

    #if (LAL_USE_OMP_SIMD == 1)
    #pragma omp simd
    #endif
    for (int i = ifrom3; i < ito3; i++) {
      v[i] += _dtfm[i] * f[i];
      x[i] += dtv * v[i];
    }

    // update angular momentum by 1/2 step
    if (igroup == 0) {
      #if (LAL_USE_OMP_SIMD == 1)
        // Workaround for compiler bug
        #ifdef __INTEL_COMPILER
        #pragma simd
        #else
        #pragma omp simd
        #endif
      #endif
      for (int i = ifrom; i < ito; i++) {
        double *quat = bonus[ellipsoid[i]].quat;
        ME_omega_richardson(dtf, dtq, angmom[i], quat, torque[i], _inertia0[i],
                            _inertia1[i], _inertia2[i]);
      }
    } else {
      #if (LAL_USE_OMP_SIMD == 1)
        // Workaround for compiler bug
        #ifdef __INTEL_COMPILER
        #pragma simd
        #else
        #pragma omp simd
        #endif
      #endif
      for (int i = ifrom; i < ito; i++) {
        if (mask[i] & groupbit) {
          double *quat = bonus[ellipsoid[i]].quat;
          ME_omega_richardson(dtf, dtq, angmom[i], quat, torque[i],
                              _inertia0[i], _inertia1[i], _inertia2[i]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereGPU::final_integrate()
{
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  double * _noalias const angmom = atom->angmom[0];
  const double * _noalias const torque = atom->torque[0];

  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (neighbor->ago == 0) {
    if (nlocal > _nlocal_max) {
      if (_nlocal_max) {
        memory->destroy(_dtfm);
        memory->destroy(_inertia0);
        memory->destroy(_inertia1);
        memory->destroy(_inertia2);
      }
      _nlocal_max = static_cast<int>(1.20 * nlocal);
      memory->create(_dtfm, _nlocal_max * 3, "fix_nve_gpu:dtfm");
      memory->create(_inertia0, _nlocal_max * 3, "fix_nve_gpu:inertia0");
      memory->create(_inertia1, _nlocal_max * 3, "fix_nve_gpu:inertia1");
      memory->create(_inertia2, _nlocal_max * 3, "fix_nve_gpu:inertia2");
    }
  }

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = comm->nthreads;
    const int tid = omp_get_thread_num();
    const int idelta = nlocal / nthreads + 1;
    const int ifrom = tid * idelta;
    const int ito = MIN(ifrom + idelta, nlocal);
    const int ifrom3 = ifrom * 3;
    const int ito3 = ito * 3;
    #else
    const int tid = 0;
    const int ifrom = 0;
    const int ifrom3 = 0;
    const int ito = nlocal;
    const int ito3 = nlocal * 3;
    #endif

    double dtfo;
    if (neighbor->ago == 0) dtfo = reset_dt_omp(ifrom, ito, tid);
    else dtfo = dtf;

    #if (LAL_USE_OMP_SIMD == 1)
    #pragma omp simd
    #endif
    for (int i = ifrom3; i < ito3; i++) {
      v[i] += _dtfm[i] * f[i];
      angmom[i] += dtfo * torque[i];
    }
  }
}

void FixNVEAsphereGPU::reset_dt() {
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (nlocal > _nlocal_max) {
    if (_nlocal_max) {
      memory->destroy(_dtfm);
      memory->destroy(_inertia0);
      memory->destroy(_inertia1);
      memory->destroy(_inertia2);
    }
    _nlocal_max = static_cast<int>(1.20 * nlocal);
    memory->create(_dtfm, _nlocal_max * 3, "fix_nve_gpu:dtfm");
    memory->create(_inertia0, _nlocal_max * 3, "fix_nve_gpu:inertia0");
    memory->create(_inertia1, _nlocal_max * 3, "fix_nve_gpu:inertia1");
    memory->create(_inertia2, _nlocal_max * 3, "fix_nve_gpu:inertia2");
  }

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = comm->nthreads;
    const int tid = omp_get_thread_num();
    const int idelta = nlocal / nthreads + 1;
    const int ifrom = tid * idelta;
    const int ito = MIN(ifrom + idelta, nlocal);
    #else
    const int tid = 0;
    const int ifrom = 0;
    const int ito = nlocal;
    #endif
    reset_dt_omp(ifrom, ito, tid);
  }
}

double FixNVEAsphereGPU::reset_dt_omp(const int ifrom, const int ito,
                                      const int tid) {
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  const int * const mask = atom->mask;

  const double dtfo = 0.5 * update->dt * force->ftm2v;
  if (tid == 0) {
    dtv = update->dt;
    dtf = dtfo;
  }

  if (igroup == 0) {
    const double * const rmass = atom->rmass;
    int n = ifrom * 3;
    for (int i = ifrom; i < ito; i++) {
      const double dtfir = dtfo / rmass[i];
      _dtfm[n++] = dtfir;
      _dtfm[n++] = dtfir;
      _dtfm[n++] = dtfir;
      double *shape = bonus[ellipsoid[i]].shape;
      double idot = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia0[i] = idot;
      idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia1[i] = idot;
      idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia2[i] = idot;
    }
  } else {
    const double * const rmass = atom->rmass;
    int n = ifrom * 3;
    for (int i = ifrom; i < ito; i++) {
      if (mask[i] & groupbit) {
        const double dtfir = dtfo / rmass[i];
        _dtfm[n++] = dtfir;
        _dtfm[n++] = dtfir;
        _dtfm[n++] = dtfir;
        double *shape = bonus[ellipsoid[i]].shape;
        double idot = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia0[i] = idot;
        idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia1[i] = idot;
        idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia2[i] = idot;
      } else {
        _dtfm[n++] = 0.0;
        _dtfm[n++] = 0.0;
        _dtfm[n++] = 0.0;
      }
    }
  }
  return dtfo;
}

double FixNVEAsphereGPU::memory_usage()
{
  return FixNVE::memory_usage() + _nlocal_max * 12 * sizeof(double);
}

