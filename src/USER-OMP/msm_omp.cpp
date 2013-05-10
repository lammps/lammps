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
  triclinic_support = 0;
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

  if (vflag_atom)
    direct_peratom<1>(n);
  else
    direct_peratom<0>(n);
}

template <int EFLAG_GLOBAL, int VFLAG_GLOBAL, int VFLAG_ATOM>
void MSMOMP::direct_eval(const int nn)
{
  const double * _noalias const * _noalias const * _noalias const qgridn = qgrid[nn];
  const double * _noalias const g_directn = g_direct[nn];
  const double * _noalias const v0_directn = v0_direct[nn];
  const double * _noalias const v1_directn = v1_direct[nn];
  const double * _noalias const v2_directn = v2_direct[nn];
  const double * _noalias const v3_directn = v3_direct[nn];
  const double * _noalias const v4_directn = v4_direct[nn];
  const double * _noalias const v5_directn = v5_direct[nn];

  double v0,v1,v2,v3,v4,v5,emsm;
  v0 = v1 = v2 = v3 = v4 = v5 = emsm = 0.0;
  const int alphan = alpha[nn];
  const int betaxn = betax[nn];
  const int betayn = betay[nn];
  const int betazn = betaz[nn];

  const int nx = nxhi_direct - nxlo_direct + 1;
  const int ny = nyhi_direct - nylo_direct + 1;

  // merge three outer loops into one for better threading

  const int nzlo_inn = nzlo_in[nn];
  const int nylo_inn = nylo_in[nn];
  const int nxlo_inn = nxlo_in[nn];
  const int numz = nzhi_in[nn] - nzlo_inn + 1;
  const int numy = nyhi_in[nn] - nylo_inn + 1;
  const int numx = nxhi_in[nn] - nxlo_inn + 1;
  const int inum = numz*numy*numx;

  const int zper = domain->zperiodic;
  const int yper = domain->yperiodic;
  const int xper = domain->xperiodic;

  const int n=nn;

#if defined(_OPENMP)
#pragma omp parallel default(none) reduction(+:v0,v1,v2,v3,v4,v5,emsm)
#endif
  {
    double esum,v0sum,v1sum,v2sum,v3sum,v4sum,v5sum;
    int i,ifrom,ito,tid,icx,icy,icz,ix,iy,iz,k;

    loop_setup_thr(ifrom, ito, tid, inum, comm->nthreads);

    for (i = ifrom; i < ito; ++i) {

      // infer outer loop indices icx, icy, icz from master loop index i

      icz = i/(numy*numx);
      icy = (i - icz*numy*numx) / numx;
      icx = i - icz*numy*numx - icy*numx;
      icz += nzlo_inn;
      icy += nylo_inn;
      icx += nxlo_inn;
      
      const int kmax = zper ? nzhi_direct : MIN(nzhi_direct,betazn - icz);
      const int jmin = yper ? nylo_direct : MAX(nylo_direct,alphan - icy);
      const int jmax = yper ? nyhi_direct : MIN(nyhi_direct,betayn - icy);
      const int imin = xper ? nxlo_direct : MAX(nxlo_direct,alphan - icx);
      const int imax = xper ? nxhi_direct : MIN(nxhi_direct,betaxn - icx);

      const double qtmp = qgridn[icz][icy][icx]; // charge on center grid point

      esum = 0.0;
      if (VFLAG_GLOBAL || VFLAG_ATOM)
        v0sum = v1sum = v2sum = v3sum = v4sum = v5sum = 0.0;

      // use hemisphere to avoid double computation of pair-wise
      //   interactions in direct sum (no computations in -z direction)

      for (iz = 1; iz <= kmax; iz++) {
        const int kk = icz+iz;
        const int zk = (iz + nzhi_direct)*ny;
        for (iy = jmin; iy <= jmax; iy++) {
          const int jj = icy+iy;
          const int zyk = (zk + iy + nyhi_direct)*nx;
          const double * _noalias const qgridnkj = &qgridn[kk][jj][icx];
          for (ix = imin; ix <= imax; ix++) {
            const double qtmp2 = qgridnkj[ix];
            k = zyk + ix + nxhi_direct;
            const double gtmp = g_directn[k];
            esum += gtmp * qtmp2;

            if (VFLAG_GLOBAL || VFLAG_ATOM) {
              v0sum += v0_directn[k] * qtmp2;
              v1sum += v1_directn[k] * qtmp2;
              v2sum += v2_directn[k] * qtmp2;
              v3sum += v3_directn[k] * qtmp2;
              v4sum += v4_directn[k] * qtmp2;
              v5sum += v5_directn[k] * qtmp2;
            }
          }
        }
      }

      // iz=0

      const int zk = nzhi_direct*ny;
      for (iy = 1; iy <= jmax; iy++) {
        const int jj = icy+iy;
        const int zyk = (zk + iy + nyhi_direct)*nx;
        const double * _noalias const qgridnkj = &qgridn[icz][jj][icx];
        for (ix = imin; ix <= imax; ix++) {
          const double qtmp2 = qgridnkj[ix];
          k = zyk + ix + nxhi_direct;
          const double gtmp = g_directn[k];
          esum += gtmp * qtmp2;

          if (VFLAG_GLOBAL || VFLAG_ATOM) {
            v0sum += v0_directn[k] * qtmp2;
            v1sum += v1_directn[k] * qtmp2;
            v2sum += v2_directn[k] * qtmp2;
            v3sum += v3_directn[k] * qtmp2;
            v4sum += v4_directn[k] * qtmp2;
            v5sum += v5_directn[k] * qtmp2;
          }
        }
      }

      // iz=0, iy=0

      const int zyk = (zk + nyhi_direct)*nx;
      const double * _noalias const qgridnkj = &qgridn[icz][icy][icx];
      for (ix = 1; ix <= imax; ix++) {
        const double qtmp2 = qgridnkj[ix];
        k = zyk + ix + nxhi_direct;
        const double gtmp = g_directn[k];
        esum += gtmp * qtmp2;

        if (VFLAG_GLOBAL || VFLAG_ATOM) {
          v0sum += v0_directn[k] * qtmp2;
          v1sum += v1_directn[k] * qtmp2;
          v2sum += v2_directn[k] * qtmp2;
          v3sum += v3_directn[k] * qtmp2;
          v4sum += v4_directn[k] * qtmp2;
          v5sum += v5_directn[k] * qtmp2;
        }
      }

      // iz=0, iy=0, ix=0

      const double qtmp2 = qgridnkj[0];
      k = zyk + nxhi_direct;
      const double gtmp = g_directn[k];
      esum += 0.5 * gtmp * qtmp2;

      // virial is zero for iz=0, iy=0, ix=0

      // accumulate per-atom energy/virial

      egrid[n][icz][icy][icx] = esum;

      if (VFLAG_ATOM) {
        v0grid[n][icz][icy][icx] = v0sum;
        v1grid[n][icz][icy][icx] = v1sum;
        v2grid[n][icz][icy][icx] = v2sum;
        v3grid[n][icz][icy][icx] = v3sum;
        v4grid[n][icz][icy][icx] = v4sum;
        v5grid[n][icz][icy][icx] = v5sum;
      }

      if (EFLAG_GLOBAL || VFLAG_GLOBAL) {
        const double qtmp3 = qgridn[icz][icy][icx];
        if (EFLAG_GLOBAL) emsm += 2.0 * esum * qtmp3;
        if (VFLAG_GLOBAL) {
          v0 += 2.0 * v0sum * qtmp3;
          v1 += 2.0 * v1sum * qtmp3;
          v2 += 2.0 * v2sum * qtmp3;
          v3 += 2.0 * v3sum * qtmp3;
          v4 += 2.0 * v4sum * qtmp3;
          v5 += 2.0 * v5sum * qtmp3;
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

template <int VFLAG_ATOM>
void MSMOMP::direct_peratom(const int nn)
{
  double * _noalias const * _noalias const * _noalias const egridn  = egrid[nn];
  double * _noalias const * _noalias const * _noalias const v0gridn  = v0grid[nn];
  double * _noalias const * _noalias const * _noalias const v1gridn  = v1grid[nn];
  double * _noalias const * _noalias const * _noalias const v2gridn  = v2grid[nn];
  double * _noalias const * _noalias const * _noalias const v3gridn  = v3grid[nn];
  double * _noalias const * _noalias const * _noalias const v4gridn  = v4grid[nn];
  double * _noalias const * _noalias const * _noalias const v5gridn  = v5grid[nn];
  const double * _noalias const * _noalias const * _noalias const qgridn = qgrid[nn];
  const double * _noalias const g_directn = g_direct[nn];
  const double * _noalias const v0_directn = v0_direct[nn];
  const double * _noalias const v1_directn = v1_direct[nn];
  const double * _noalias const v2_directn = v2_direct[nn];
  const double * _noalias const v3_directn = v3_direct[nn];
  const double * _noalias const v4_directn = v4_direct[nn];
  const double * _noalias const v5_directn = v5_direct[nn];


  const int alphan = alpha[nn];
  const int betaxn = betax[nn];
  const int betayn = betay[nn];
  const int betazn = betaz[nn];

  const int nx = nxhi_direct - nxlo_direct + 1;
  const int ny = nyhi_direct - nylo_direct + 1;

  // merge three outer loops into one

  const int nzlo_inn = nzlo_in[nn];
  const int nylo_inn = nylo_in[nn];
  const int nxlo_inn = nxlo_in[nn];
  const int numz = nzhi_in[nn] - nzlo_inn + 1;
  const int numy = nyhi_in[nn] - nylo_inn + 1;
  const int numx = nxhi_in[nn] - nxlo_inn + 1;
  const int inum = numz*numy*numx;

  const int zper = domain->zperiodic;
  const int yper = domain->yperiodic;
  const int xper = domain->xperiodic;

  const int n=nn;
  int i,ifrom,ito,tid,icx,icy,icz,ix,iy,iz,k;


  for (i = 0; i < inum; ++i) {

    // infer outer loop indices icx, icy, icz from master loop index i

    icz = i/(numy*numx);
    icy = (i - icz*numy*numx) / numx;
    icx = i - icz*numy*numx - icy*numx;
    icz += nzlo_inn;
    icy += nylo_inn;
    icx += nxlo_inn;
    
    const int kmax = zper ? nzhi_direct : MIN(nzhi_direct,betazn - icz);
    const int jmin = yper ? nylo_direct : MAX(nylo_direct,alphan - icy);
    const int jmax = yper ? nyhi_direct : MIN(nyhi_direct,betayn - icy);
    const int imin = xper ? nxlo_direct : MAX(nxlo_direct,alphan - icx);
    const int imax = xper ? nxhi_direct : MIN(nxhi_direct,betaxn - icx);

    const double qtmp = qgridn[icz][icy][icx]; // charge on center grid point


    // use hemisphere to avoid double computation of pair-wise
    //   interactions in direct sum (no computations in -z direction)

    for (iz = 1; iz <= kmax; iz++) {
      const int kk = icz+iz;
      const int zk = (iz + nzhi_direct)*ny;
      for (iy = jmin; iy <= jmax; iy++) {
        const int jj = icy+iy;
        const int zyk = (zk + iy + nyhi_direct)*nx;
        double * _noalias const egridnkj = &egridn[kk][jj][icx];
        for (ix = imin; ix <= imax; ix++) {
          k = zyk + ix + nxhi_direct;
          const int ii = icx+ix;
          const double gtmp = g_directn[k];

          egridnkj[ix] += gtmp * qtmp;

          if (VFLAG_ATOM) {
            v0gridn[kk][jj][ii] += v0_directn[k] * qtmp;
            v1gridn[kk][jj][ii] += v1_directn[k] * qtmp;
            v2gridn[kk][jj][ii] += v2_directn[k] * qtmp;
            v3gridn[kk][jj][ii] += v3_directn[k] * qtmp;
            v4gridn[kk][jj][ii] += v4_directn[k] * qtmp;
            v5gridn[kk][jj][ii] += v5_directn[k] * qtmp;
          }
        }
      }
    }

    // iz=0

    const int zk = nzhi_direct*ny;
    for (iy = 1; iy <= jmax; iy++) {
      const int jj = icy+iy;
      const int zyk = (zk + iy + nyhi_direct)*nx;
      double * _noalias const egridnkj = &egridn[icz][jj][icx];
      for (ix = imin; ix <= imax; ix++) {
        k = zyk + ix + nxhi_direct;
        const int ii = icx+ix;
        const double gtmp = g_directn[k];

        egridnkj[ix] += gtmp * qtmp;

        if (VFLAG_ATOM) {
          v0gridn[icz][jj][ii] += v0_directn[k] * qtmp;
          v1gridn[icz][jj][ii] += v1_directn[k] * qtmp;
          v2gridn[icz][jj][ii] += v2_directn[k] * qtmp;
          v3gridn[icz][jj][ii] += v3_directn[k] * qtmp;
          v4gridn[icz][jj][ii] += v4_directn[k] * qtmp;
          v5gridn[icz][jj][ii] += v5_directn[k] * qtmp;
        }
      }
    }

    // iz=0, iy=0

    const int zyk = (zk + nyhi_direct)*nx;
    double * _noalias const egridnkj = &egridn[icz][icy][icx];
    for (ix = 1; ix <= imax; ix++) {
      k = zyk + ix + nxhi_direct;
      const int ii = icx+ix;
      const double gtmp = g_directn[k];

      egridnkj[ix] += gtmp * qtmp;

      if (VFLAG_ATOM) {
        v0gridn[icz][icy][ii] += v0_directn[k] * qtmp;
        v1gridn[icz][icy][ii] += v1_directn[k] * qtmp;
        v2gridn[icz][icy][ii] += v2_directn[k] * qtmp;
        v3gridn[icz][icy][ii] += v3_directn[k] * qtmp;
        v4gridn[icz][icy][ii] += v4_directn[k] * qtmp;
        v5gridn[icz][icy][ii] += v5_directn[k] * qtmp;
      }
    }

    // iz=0, iy=0, ix=0

    k = zyk + nxhi_direct;
    const double gtmp = g_directn[k];
    egridnkj[0] += 0.5 * gtmp * qtmp;

    // virial is zero for iz=0, iy=0, ix=0

  }
}
