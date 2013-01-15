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
   Contributing authors: Axel Kohlmeyer (Temple U), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "msm_omp.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "math_const.h"

#include <string.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "suffix.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MSMOMP::MSMOMP(LAMMPS *lmp, int narg, char **arg) :
  MSM(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE)
{
  suffix_flag |= Suffix::OMP;
}

/* ----------------------------------------------------------------------
   run the regular toplevel compute method from plain PPPPM
   which will have individual methods replaced by our threaded
   versions and then call the obligatory force reduction.
------------------------------------------------------------------------- */

void MSMOMP::compute(int eflag, int vflag)
{

  MSM::compute(eflag,vflag);

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(eflag,vflag)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ----------------------------------------------------------------------
   MSM direct part procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSMOMP::direct(int n) 
{
  // zero out electric potential

  memset(&(egrid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
  
  // zero out virial

  if (vflag_atom) {
    memset(&(v0grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v1grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v2grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v3grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v4grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
    memset(&(v5grid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));
  }

  if (eflag_global) {
    if (vflag_global) {
      if (vflag_atom)
        direct_eval<1,1,1>(n);
      else
        direct_eval<1,1,0>(n);
    } else {
      if (vflag_atom)
        direct_eval<1,0,1>(n);
      else
        direct_eval<1,0,0>(n);
    }
  } else { // !eflag_global
    if (vflag_global) {
      if (vflag_atom)
        direct_eval<0,1,1>(n);
      else
        direct_eval<0,1,0>(n);
    } else {
      if (vflag_atom)
        direct_eval<0,0,1>(n);
      else
        direct_eval<0,0,0>(n);
    }
  }
}

template <int EFLAG_GLOBAL, int VFLAG_GLOBAL, int VFLAG_ATOM>
void MSMOMP::direct_eval(const int n)
{
  double * const * const * const egridn  = egrid[n];
  const double * const * const * const qgridn = qgrid[n];
  const double * const g_directn = g_direct[n];
  const double * const v0_directn = v0_direct[n];
  const double * const v1_directn = v1_direct[n];
  const double * const v2_directn = v2_direct[n];
  const double * const v3_directn = v3_direct[n];
  const double * const v4_directn = v4_direct[n];
  const double * const v5_directn = v5_direct[n];

  double v0,v1,v2,v3,v4,v5,emsm;
  v0 = v1 = v2 = v3 = v4 = v5 = emsm = 0.0;
  const int alphan = alpha[n];
  const int betaxn = betax[n];
  const int betayn = betay[n];
  const int betazn = betaz[n];

  const int nx = nxhi_direct - nxlo_direct + 1;
  const int ny = nyhi_direct - nylo_direct + 1;

  // merge three outer loops into one for better threading

  const int nzlo_inn = nzlo_in[n];
  const int nylo_inn = nylo_in[n];
  const int nxlo_inn = nxlo_in[n];
  const int numz = nzhi_in[n] - nzlo_inn + 1;
  const int numy = nyhi_in[n] - nylo_inn + 1;
  const int numx = nxhi_in[n] - nxlo_inn + 1;
  const int inum = numz*numy*numx;

  const int zper = domain->zperiodic;
  const int yper = domain->yperiodic;
  const int xper = domain->xperiodic;

#if defined(_OPENMP)
#pragma omp parallel default(none) reduction(+:v0,v1,v2,v3,v4,v5,emsm)
#endif
  {
    double qtmp,esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
    int i,ifrom,ito,tid,icx,icy,icz,ix,iy,iz,k;

    loop_setup_thr(ifrom, ito, tid, inum, comm->nthreads);

    for (i = ifrom; i < ito; ++i) {

      // infer outer loop indices icx, icy, icz from master loop index

      icz = i/(numy*numx);
      icy = (i - icz*numy*numx) / numx;
      icx = i - icz*numy*numx - icy*numx;
      icz += nzlo_inn;
      icy += nylo_inn;
      icx += nxlo_inn;
      
      const int kmin = zper ? nzlo_direct : MAX(nzlo_direct,alphan - icz);
      const int kmax = zper ? nzhi_direct : MIN(nzhi_direct,betazn - icz);
      const int jmin = yper ? nylo_direct : MAX(nylo_direct,alphan - icy);
      const int jmax = yper ? nyhi_direct : MIN(nyhi_direct,betayn - icy);
      const int imin = xper ? nxlo_direct : MAX(nxlo_direct,alphan - icx);
      const int imax = xper ? nxhi_direct : MIN(nxhi_direct,betaxn - icx);

      esum = 0.0;
      if (VFLAG_GLOBAL || VFLAG_ATOM)
        v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

      for (iz = kmin; iz <= kmax; iz++) {
        const int kk = icz+iz;
        const int zk = (iz + nzhi_direct)*ny;
        for (iy = jmin; iy <= jmax; iy++) {
          const int jj = icy+iy;
          const int zyk = (zk + iy + nyhi_direct)*nx;
          for (ix = imin; ix <= imax; ix++) {
            qtmp = qgridn[kk][jj][icx+ix];
            k = zyk + ix + nxhi_direct;
            esum += g_directn[k] * qtmp;

            if (VFLAG_GLOBAL || VFLAG_ATOM) {
              v0sum += v0_directn[k] * qtmp;
              v1sum += v1_directn[k] * qtmp;
              v2sum += v2_directn[k] * qtmp;
              v3sum += v3_directn[k] * qtmp;
              v4sum += v4_directn[k] * qtmp;
              v5sum += v5_directn[k] * qtmp;
            }
          }
        }
      }
      egridn[icz][icy][icx] = esum;

      if (VFLAG_ATOM) {
        v0grid[n][icz][icy][icx] = v0sum;
        v1grid[n][icz][icy][icx] = v1sum;
        v2grid[n][icz][icy][icx] = v2sum;
        v3grid[n][icz][icy][icx] = v3sum;
        v4grid[n][icz][icy][icx] = v4sum;
        v5grid[n][icz][icy][icx] = v5sum;
      }

      if (EFLAG_GLOBAL || VFLAG_GLOBAL) {
        qtmp = qgridn[icz][icy][icx];
        if (EFLAG_GLOBAL) emsm += esum * qtmp;
        if (VFLAG_GLOBAL) {
          v0 += v0sum * qtmp;
          v1 += v1sum * qtmp;
          v2 += v2sum * qtmp;
          v3 += v3sum * qtmp;
          v4 += v4sum * qtmp;
          v5 += v5sum * qtmp;
        }
      }
    }
  } // end of omp parallel region

  if (EFLAG_GLOBAL || VFLAG_GLOBAL) {
    if (EFLAG_GLOBAL) energy += emsm;
    if (VFLAG_GLOBAL) {
      virial[0] += v0;
      virial[1] += v1;
      virial[2] += v2;
      virial[3] += v3;
      virial[4] += v4;
      virial[5] += v5;
    }
  }
}
