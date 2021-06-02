// clang-format off
/* -------------------------------------------------------------------------
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
   Contributing author: Axel Kohlmeyer (Temple U)
   per-thread data management for LAMMPS
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdio>

#include "thr_data.h"

#include "memory.h"
#include "timer.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ThrData::ThrData(int tid, Timer *t)
  : _f(0),_torque(0),_erforce(0),_de(0),_drho(0),_mu(0),_lambda(0),_rhoB(0),
    _D_values(0),_rho(0),_fp(0),_rho1d(0),_drho1d(0),_rho1d_6(0),_drho1d_6(0),
    _tid(tid), _timer(t)
{
  _timer_active = 0;
}


/* ---------------------------------------------------------------------- */

void ThrData::check_tid(int tid)
{
  if (tid != _tid)
    fprintf(stderr,"WARNING: external and internal tid mismatch %d != %d\n",tid,_tid);
}

/* ---------------------------------------------------------------------- */

void ThrData::_stamp(enum Timer::ttype flag)
{
  // do nothing until it gets set to 0 in ::setup()
  if (_timer_active < 0) return;

  if (flag == Timer::START) {
    _timer_active = 1;
  }

  if (_timer_active) _timer->stamp(flag);
}

/* ---------------------------------------------------------------------- */

double ThrData::get_time(enum Timer::ttype flag)
{
  if (_timer)
    return _timer->get_wall(flag);
  else
    return 0.0;
}

/* ---------------------------------------------------------------------- */

void ThrData::init_force(int nall, double **f, double **torque,
                         double *erforce, double *de, double *drho)
{
  eng_vdwl=eng_coul=eng_bond=eng_angle=eng_dihed=eng_imprp=eng_kspce=0.0;
  memset(virial_pair,0,6*sizeof(double));
  memset(virial_bond,0,6*sizeof(double));
  memset(virial_angle,0,6*sizeof(double));
  memset(virial_dihed,0,6*sizeof(double));
  memset(virial_imprp,0,6*sizeof(double));
  memset(virial_kspce,0,6*sizeof(double));

  eatom_pair=eatom_bond=eatom_angle=eatom_dihed=eatom_imprp=eatom_kspce=nullptr;
  vatom_pair=vatom_bond=vatom_angle=vatom_dihed=vatom_imprp=vatom_kspce=nullptr;

  if (nall >= 0 && f) {
    _f = f + _tid*nall;
    memset(&(_f[0][0]),0,nall*3*sizeof(double));
  } else _f = nullptr;

  if (nall >= 0 && torque) {
    _torque = torque + _tid*nall;
    memset(&(_torque[0][0]),0,nall*3*sizeof(double));
  } else _torque = nullptr;

  if (nall >= 0 && erforce) {
    _erforce = erforce + _tid*nall;
    memset(&(_erforce[0]),0,nall*sizeof(double));
  } else _erforce = nullptr;

  if (nall >= 0 && de) {
    _de = de + _tid*nall;
    memset(&(_de[0]),0,nall*sizeof(double));
  } else _de = nullptr;

  if (nall >= 0 && drho) {
    _drho = drho + _tid*nall;
    memset(&(_drho[0]),0,nall*sizeof(double));
  } else _drho = nullptr;
}

/* ----------------------------------------------------------------------
   set up and clear out locally managed per atom arrays
------------------------------------------------------------------------- */

void ThrData::init_eam(int nall, double *rho)
{
  if (nall >= 0 && rho) {
    _rho = rho + _tid*nall;
    memset(_rho, 0, nall*sizeof(double));
  }
}

/* ---------------------------------------------------------------------- */

void ThrData::init_adp(int nall, double *rho, double **mu, double **lambda)
{
  init_eam(nall, rho);

  if (nall >= 0 && mu && lambda) {
    _mu = mu + _tid*nall;
    _lambda = lambda + _tid*nall;
    memset(&(_mu[0][0]), 0, nall*3*sizeof(double));
    memset(&(_lambda[0][0]), 0, nall*6*sizeof(double));
  }
}

/* ---------------------------------------------------------------------- */

void ThrData::init_eim(int nall, double *rho, double *fp)
{
  init_eam(nall, rho);

  if (nall >= 0 && fp) {
    _fp = fp + _tid*nall;
    memset(_fp,0,nall*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   if order > 0 : set up per thread storage for PPPM
   if order < 0 : free per thread storage for PPPM
------------------------------------------------------------------------- */
#if defined(FFT_SINGLE)
typedef float FFT_SCALAR;
#else
typedef double FFT_SCALAR;
#endif

void ThrData::init_pppm(int order, Memory *memory)
{
  FFT_SCALAR **rho1d, **drho1d;
  if (order > 0) {
    rho1d = static_cast<FFT_SCALAR **>(_rho1d);
    drho1d = static_cast<FFT_SCALAR **>(_drho1d);
    if (rho1d) memory->destroy2d_offset(rho1d,-order/2);
    if (drho1d) memory->destroy2d_offset(drho1d,-order/2);
    memory->create2d_offset(rho1d,3,-order/2,order/2,"thr_data:rho1d");
    memory->create2d_offset(drho1d,3,-order/2,order/2,"thr_data:drho1d");
    _rho1d = static_cast<void *>(rho1d);
    _drho1d = static_cast<void *>(drho1d);
  } else {
    order = -order;
    rho1d = static_cast<FFT_SCALAR **>(_rho1d);
    drho1d = static_cast<FFT_SCALAR **>(_drho1d);
    if (rho1d) memory->destroy2d_offset(rho1d,-order/2);
    if (drho1d) memory->destroy2d_offset(drho1d,-order/2);
    _rho1d = nullptr;
    _drho1d = nullptr;
  }
}

/* ----------------------------------------------------------------------
   if order > 0 : set up per thread storage for PPPM
   if order < 0 : free per thread storage for PPPM
------------------------------------------------------------------------- */
#if defined(FFT_SINGLE)
typedef float FFT_SCALAR;
#else
typedef double FFT_SCALAR;
#endif

void ThrData::init_pppm_disp(int order_6, Memory *memory)
{
  FFT_SCALAR **rho1d_6, **drho1d_6;
  if (order_6 > 0) {
    rho1d_6 = static_cast<FFT_SCALAR **>(_rho1d_6);
    drho1d_6 = static_cast<FFT_SCALAR **>(_drho1d_6);
    if (rho1d_6) memory->destroy2d_offset(rho1d_6,-order_6/2);
    if (drho1d_6) memory->destroy2d_offset(drho1d_6,-order_6/2);
    memory->create2d_offset(rho1d_6,3,-order_6/2,order_6/2,"thr_data:rho1d_6");
    memory->create2d_offset(drho1d_6,3,-order_6/2,order_6/2,"thr_data:drho1d_6");
    _rho1d_6 = static_cast<void *>(rho1d_6);
    _drho1d_6 = static_cast<void *>(drho1d_6);
  } else {
    order_6 = -order_6;
    rho1d_6 = static_cast<FFT_SCALAR **>(_rho1d_6);
    drho1d_6 = static_cast<FFT_SCALAR **>(_drho1d_6);
    if (rho1d_6) memory->destroy2d_offset(rho1d_6,-order_6/2);
    if (drho1d_6) memory->destroy2d_offset(drho1d_6,-order_6/2);
  }
}

/* ----------------------------------------------------------------------
   compute global pair virial via summing F dot r over own & ghost atoms
   at this point, only pairwise forces have been accumulated in atom->f
------------------------------------------------------------------------- */

void ThrData::virial_fdotr_compute(double **x, int nlocal, int nghost, int nfirst)
{

  // sum over force on all particles including ghosts

  if (nfirst < 0) {
    int nall = nlocal + nghost;
    for (int i = 0; i < nall; i++) {
      virial_pair[0] += _f[i][0]*x[i][0];
      virial_pair[1] += _f[i][1]*x[i][1];
      virial_pair[2] += _f[i][2]*x[i][2];
      virial_pair[3] += _f[i][1]*x[i][0];
      virial_pair[4] += _f[i][2]*x[i][0];
      virial_pair[5] += _f[i][2]*x[i][1];
    }

  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = nfirst;
    for (int i = 0; i < nall; i++) {
      virial_pair[0] += _f[i][0]*x[i][0];
      virial_pair[1] += _f[i][1]*x[i][1];
      virial_pair[2] += _f[i][2]*x[i][2];
      virial_pair[3] += _f[i][1]*x[i][0];
      virial_pair[4] += _f[i][2]*x[i][0];
      virial_pair[5] += _f[i][2]*x[i][1];
    }
    nall = nlocal + nghost;
    for (int i = nlocal; i < nall; i++) {
      virial_pair[0] += _f[i][0]*x[i][0];
      virial_pair[1] += _f[i][1]*x[i][1];
      virial_pair[2] += _f[i][2]*x[i][2];
      virial_pair[3] += _f[i][1]*x[i][0];
      virial_pair[4] += _f[i][2]*x[i][0];
      virial_pair[5] += _f[i][2]*x[i][1];
    }
  }
}

/* ---------------------------------------------------------------------- */

double ThrData::memory_usage()
{
  double bytes = (7 + 6*6) * sizeof(double);
  bytes += (double)2 * sizeof(double*);
  bytes += (double)4 * sizeof(int);

  return bytes;
}

/* additional helper functions */

// reduce per thread data into the first part of the data
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on
// multi-dimensional arrays like force stored in this order
// x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with writing to the array .
void LAMMPS_NS::data_reduce_thr(double *dall, int nall, int nthreads, int ndim, int tid)
{
#if defined(_OPENMP)
  // NOOP in single-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    const int nvals = ndim*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nvals) ? nvals : (ifrom + idelta);

#if defined(USER_OMP_NO_UNROLL)
    if (ifrom < nvals) {
      int m = 0;

      for (m = ifrom; m < ito; ++m) {
        for (int n = 1; n < nthreads; ++n) {
          dall[m] += dall[n*nvals + m];
          dall[n*nvals + m] = 0.0;
        }
      }
    }
#else
    // this if protects against having more threads than atoms
    if (ifrom < nvals) {
      int m = 0;

      // for architectures that have L1 D-cache line sizes of 64 bytes
      // (8 doubles) wide, explicitly unroll this loop to  compute 8
      // contiguous values in the array at a time
      // -- modify this code based on the size of the cache line
      double t0, t1, t2, t3, t4, t5, t6, t7;
      for (m = ifrom; m < (ito-7); m+=8) {
        t0 = dall[m  ];
        t1 = dall[m+1];
        t2 = dall[m+2];
        t3 = dall[m+3];
        t4 = dall[m+4];
        t5 = dall[m+5];
        t6 = dall[m+6];
        t7 = dall[m+7];
        for (int n = 1; n < nthreads; ++n) {
          t0 += dall[n*nvals + m  ];
          t1 += dall[n*nvals + m+1];
          t2 += dall[n*nvals + m+2];
          t3 += dall[n*nvals + m+3];
          t4 += dall[n*nvals + m+4];
          t5 += dall[n*nvals + m+5];
          t6 += dall[n*nvals + m+6];
          t7 += dall[n*nvals + m+7];
          dall[n*nvals + m  ] = 0.0;
          dall[n*nvals + m+1] = 0.0;
          dall[n*nvals + m+2] = 0.0;
          dall[n*nvals + m+3] = 0.0;
          dall[n*nvals + m+4] = 0.0;
          dall[n*nvals + m+5] = 0.0;
          dall[n*nvals + m+6] = 0.0;
          dall[n*nvals + m+7] = 0.0;
        }
        dall[m  ] = t0;
        dall[m+1] = t1;
        dall[m+2] = t2;
        dall[m+3] = t3;
        dall[m+4] = t4;
        dall[m+5] = t5;
        dall[m+6] = t6;
        dall[m+7] = t7;
      }
      // do the last < 8 values
      for (; m < ito; m++) {
        for (int n = 1; n < nthreads; ++n) {
          dall[m] += dall[n*nvals + m];
          dall[n*nvals + m] = 0.0;
        }
      }
    }
#endif
  }
#else
  // NOOP in non-threaded execution.
  return;
#endif
}
