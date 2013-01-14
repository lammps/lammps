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
   MSM direct part procedure for intermediate grid levels
------------------------------------------------------------------------- */

void MSMOMP::direct(int n) 
{
  // zero out electric potential

  memset(&(egrid[n][nzlo_out[n]][nylo_out[n]][nxlo_out[n]]),0,ngrid[n]*sizeof(double));

  if (evflag) direct_eval<1>(n);
  else direct_eval<0>(n);
}

template <int EVFLAG> void MSMOMP::direct_eval(const int n)
{
  //fprintf(screen,"Direct contribution on level %i\n\n",n);

  double *** egridn = egrid[n];
  const double * const * const * const qgridn = qgrid[n];
  double v0,v1,v2,v3,v4,v5,emsm;
  v0 = v1 = v2 = v3 = v4 = v5 = emsm = 0.0;
  const int alphan = alpha[n];
  const int betaxn = betax[n];
  const int betayn = betay[n];
  const int betazn = betaz[n];

#if defined(_OPENMP)
#pragma omp parallel default(none) shared(egridn) reduction(+:v0,v1,v2,v3,v4,v5,emsm)
#endif
  {
    double qtmp,etmp,esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
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

          if (EVFLAG)
            esum = v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

          etmp = 0.0;
          for (iz = kmin; iz <= kmax; iz++) {
            kk = icz+iz;
            zk = (iz + nzhi_direct)*ny;
            for (iy = jmin; iy <= jmax; iy++) {
              jj = icy+iy;
              zyk = (zk + iy + nyhi_direct)*nx;
              for (ix = imin; ix <= imax; ix++) {
                qtmp = qgridn[kk][jj][icx+ix];
                k = zyk + ix + nxhi_direct;
                etmp += g_direct[n][k] * qtmp;

                if (EVFLAG) {
                  if (eflag_global) esum += g_direct[n][k] * qtmp;
                  if (vflag_global) {
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
          }
          egridn[icz][icy][icx] = etmp;


          if (EVFLAG) {
            qtmp = qgridn[icz][icy][icx];
            if (eflag_global) emsm += esum * qtmp;
            if (vflag_global) {
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
  }
  if (EVFLAG) {
    if (eflag_global) energy += emsm;
    if (vflag_global) {
      virial[0] += v0;
      virial[1] += v1;
      virial[2] += v2;
      virial[3] += v3;
      virial[4] += v4;
      virial[5] += v5;
    }
  }
}

