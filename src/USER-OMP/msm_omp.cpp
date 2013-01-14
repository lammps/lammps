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
   Contributing author: Axel Kohlmeyer (Temple U)
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

  if (evflag) {
    if (eflag_global) {
      if (vflag_global)
        direct_eval<1,1,1>(n);
      else
        direct_eval<1,1,0>(n);
    } else {
      if (vflag_global)
        direct_eval<1,0,1>(n);
      else
        direct_eval<1,0,0>(n);
    }
  } else {
    direct_eval<0,0,0>(n);
  }
}

template <int EVFLAG, int EFLAG_GLOBAL, int VFLAG_GLOBAL>
void MSMOMP::direct_eval(const int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double * const * const * const egridn = egrid[n];
  const double * const * const * const qgridn = qgrid[n];
  const double * const g_directn = g_direct[n];

  double v0,v1,v2,v3,v4,v5,emsm;
  v0 = v1 = v2 = v3 = v4 = v5 = emsm = 0.0;
  const int alphan = alpha[n];
  const int betaxn = betax[n];
  const int betayn = betay[n];
  const int betazn = betaz[n];

#if defined(_OPENMP)
#pragma omp parallel default(none) reduction(+:v0,v1,v2,v3,v4,v5,emsm)
#endif
  {
    double qtmp,esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
    int icx,icy,icz,ix,iy,iz,zk,zyk,k;
    int jj,kk;
    int imin,imax,jmin,jmax,kmin,kmax;
  
    const int nx = nxhi_direct - nxlo_direct + 1;
    const int ny = nyhi_direct - nylo_direct + 1;

#if defined(_OPENMP)
    const int tid = omp_get_thread_num();

    // each thread works on a fixed chunk of grid planes
    const int inum = nzhi_in[n] - nzlo_in[n] + 1;
    const int idelta = 1 + inum/comm->nthreads;
    const int iczfrom = nzlo_in[n] + tid*idelta;
    const int iczto   = ((iczfrom + idelta) > nzhi_in[n]+1) ? nzhi_in[n]+1 : iczfrom + idelta;
#else
    const int tid = 0;
    const int iczfrom = nzlo_in[n];
    const int iczto = nzhi_in[n];
#endif

    for (icz = iczfrom; icz < iczto; icz++) {

      if (domain->zperiodic) {
        kmin = nzlo_direct;
        kmax = nzhi_direct;
      } else {
        kmin = MAX(nzlo_direct,alphan - icz);
        kmax = MIN(nzhi_direct,betazn - icz);
      }
        
      for (icy = nylo_in[n]; icy <= nyhi_in[n]; icy++) {

        if (domain->yperiodic) {
          jmin = nylo_direct;
          jmax = nyhi_direct;
        } else {
          jmin = MAX(nylo_direct,alphan - icy);
          jmax = MIN(nyhi_direct,betayn - icy);
        }
        
        for (icx = nxlo_in[n]; icx <= nxhi_in[n]; icx++) {

          if (domain->xperiodic) {
            imin = nxlo_direct;
            imax = nxhi_direct;
          } else {
            imin = MAX(nxlo_direct,alphan - icx);
            imax = MIN(nxhi_direct,betaxn - icx);
          }

          if (VFLAG_GLOBAL)
            v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

          esum = 0.0;
          for (iz = kmin; iz <= kmax; iz++) {
            kk = icz+iz;
            zk = (iz + nzhi_direct)*ny;
            for (iy = jmin; iy <= jmax; iy++) {
              jj = icy+iy;
              zyk = (zk + iy + nyhi_direct)*nx;
              for (ix = imin; ix <= imax; ix++) {
                qtmp = qgridn[kk][jj][icx+ix];
                k = zyk + ix + nxhi_direct;
                esum += g_directn[k] * qtmp;

                if (VFLAG_GLOBAL) {
                    v0sum += v0_direct[n][k] * qtmp;
                    v1sum += v1_direct[n][k] * qtmp;
                    v2sum += v2_direct[n][k] * qtmp;
                    v3sum += v3_direct[n][k] * qtmp;
                    v4sum += v4_direct[n][k] * qtmp;
                    v5sum += v5_direct[n][k] * qtmp;
                }
              }
            }
          }
          egridn[icz][icy][icx] = esum;

          if (EVFLAG) {
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
      }
    }
  } // end of omp parallel region

  if (EVFLAG) {
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
