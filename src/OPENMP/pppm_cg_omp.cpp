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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pppm_cg_omp.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

#define EPS_HOC 1.0e-7

/* ---------------------------------------------------------------------- */

PPPMCGOMP::PPPMCGOMP(LAMMPS *lmp) : PPPMCG(lmp), ThrOMP(lmp, THR_KSPACE)
{
  triclinic_support = 0;
  suffix_flag |= Suffix::OMP;
}

/* ----------------------------------------------------------------------
   clean up per-thread allocations
------------------------------------------------------------------------- */

PPPMCGOMP::~PPPMCGOMP()
{
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    thr->init_pppm(-order,memory);
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMCGOMP::allocate()
{
  PPPMCG::allocate();

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    thr->init_pppm(order,memory);
  }
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void PPPMCGOMP::compute_gf_ik()
{
  const double * const prd = (triclinic==0) ? domain->prd : domain->prd_lamda;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  const int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) *
                                    pow(-log(EPS_HOC),0.25));
  const int numk = nxhi_fft - nxlo_fft + 1;
  const int numl = nyhi_fft - nylo_fft + 1;

  const int twoorder = 2*order;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    double snx,sny,snz;
    double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
    double sum1,dot1,dot2;
    double numerator,denominator;
    double sqk;

    int k,l,m,nx,ny,nz,kper,lper,mper,n,nfrom,nto,tid;

    loop_setup_thr(nfrom, nto, tid, nfft, comm->nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);

    for (n = nfrom; n < nto; ++n) {
      m = n / (numl*numk);
      l = (n - m*numl*numk) / numk;
      k = n - m*numl*numk - l*numk;
      m += nzlo_fft;
      l += nylo_fft;
      k += nxlo_fft;

      mper = m - nz_pppm*(2*m/nz_pppm);
      snz = square(sin(0.5*unitkz*mper*zprd_slab/nz_pppm));

      lper = l - ny_pppm*(2*l/ny_pppm);
      sny = square(sin(0.5*unitky*lper*yprd/ny_pppm));

      kper = k - nx_pppm*(2*k/nx_pppm);
      snx = square(sin(0.5*unitkx*kper*xprd/nx_pppm));

      sqk = square(unitkx*kper) + square(unitky*lper) + square(unitkz*mper);

      if (sqk != 0.0) {
        numerator = 12.5663706/sqk;
        denominator = gf_denom(snx,sny,snz);
        sum1 = 0.0;

        for (nx = -nbx; nx <= nbx; nx++) {
          qx = unitkx*(kper+nx_pppm*nx);
          sx = exp(-0.25*square(qx/g_ewald));
          argx = 0.5*qx*xprd/nx_pppm;
          wx = powsinxx(argx,twoorder);

          for (ny = -nby; ny <= nby; ny++) {
            qy = unitky*(lper+ny_pppm*ny);
            sy = exp(-0.25*square(qy/g_ewald));
            argy = 0.5*qy*yprd/ny_pppm;
            wy = powsinxx(argy,twoorder);

            for (nz = -nbz; nz <= nbz; nz++) {
              qz = unitkz*(mper+nz_pppm*nz);
              sz = exp(-0.25*square(qz/g_ewald));
              argz = 0.5*qz*zprd_slab/nz_pppm;
              wz = powsinxx(argz,twoorder);

              dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
              dot2 = qx*qx+qy*qy+qz*qz;
              sum1 += (dot1/dot2) * sx*sy*sz * wx*wy*wz;
            }
          }
        }
        greensfn[n] = numerator*sum1/denominator;
      } else greensfn[n] = 0.0;
    }
    thr->timer(Timer::KSPACE);
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
------------------------------------------------------------------------- */

void PPPMCGOMP::compute_gf_ad()
{

  const double * const prd = (triclinic==0) ? domain->prd : domain->prd_lamda;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  const double unitkx = (MY_2PI/xprd);
  const double unitky = (MY_2PI/yprd);
  const double unitkz = (MY_2PI/zprd_slab);

  const int numk = nxhi_fft - nxlo_fft + 1;
  const int numl = nyhi_fft - nylo_fft + 1;

  const int twoorder = 2*order;
  double sf0=0.0,sf1=0.0,sf2=0.0,sf3=0.0,sf4=0.0,sf5=0.0;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE reduction(+:sf0,sf1,sf2,sf3,sf4,sf5)
#endif
  {
    double snx,sny,snz,sqk;
    double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
    double numerator,denominator;
    int k,l,m,kper,lper,mper,n,nfrom,nto,tid;

    loop_setup_thr(nfrom, nto, tid, nfft, comm->nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);

    for (n = nfrom; n < nto; ++n) {

      m = n / (numl*numk);
      l = (n - m*numl*numk) / numk;
      k = n - m*numl*numk - l*numk;
      m += nzlo_fft;
      l += nylo_fft;
      k += nxlo_fft;

      mper = m - nz_pppm*(2*m/nz_pppm);
      qz = unitkz*mper;
      snz = square(sin(0.5*qz*zprd_slab/nz_pppm));
      sz = exp(-0.25*square(qz/g_ewald));
      argz = 0.5*qz*zprd_slab/nz_pppm;
      wz = powsinxx(argz,twoorder);

      lper = l - ny_pppm*(2*l/ny_pppm);
      qy = unitky*lper;
      sny = square(sin(0.5*qy*yprd/ny_pppm));
      sy = exp(-0.25*square(qy/g_ewald));
      argy = 0.5*qy*yprd/ny_pppm;
      wy = powsinxx(argy,twoorder);

      kper = k - nx_pppm*(2*k/nx_pppm);
      qx = unitkx*kper;
      snx = square(sin(0.5*qx*xprd/nx_pppm));
      sx = exp(-0.25*square(qx/g_ewald));
      argx = 0.5*qx*xprd/nx_pppm;
      wx = powsinxx(argx,twoorder);

      sqk = qx*qx + qy*qy + qz*qz;

      if (sqk != 0.0) {
        numerator = MY_4PI/sqk;
        denominator = gf_denom(snx,sny,snz);
        greensfn[n] = numerator*sx*sy*sz*wx*wy*wz/denominator;
        sf0 += sf_precoeff1[n]*greensfn[n];
        sf1 += sf_precoeff2[n]*greensfn[n];
        sf2 += sf_precoeff3[n]*greensfn[n];
        sf3 += sf_precoeff4[n]*greensfn[n];
        sf4 += sf_precoeff5[n]*greensfn[n];
        sf5 += sf_precoeff6[n]*greensfn[n];
      } else {
        greensfn[n] = 0.0;
        sf0 += sf_precoeff1[n]*greensfn[n];
        sf1 += sf_precoeff2[n]*greensfn[n];
        sf2 += sf_precoeff3[n]*greensfn[n];
        sf3 += sf_precoeff4[n]*greensfn[n];
        sf4 += sf_precoeff5[n]*greensfn[n];
        sf5 += sf_precoeff6[n]*greensfn[n];
      }
    }
    thr->timer(Timer::KSPACE);
  } // end of parallel region

  // compute the coefficients for the self-force correction

  double prex, prey, prez, tmp[6];
  prex = prey = prez = MY_PI/volume;
  prex *= nx_pppm/xprd;
  prey *= ny_pppm/yprd;
  prez *= nz_pppm/zprd_slab;
  tmp[0] = sf0 * prex;
  tmp[1] = sf1 * prex*2;
  tmp[2] = sf2 * prey;
  tmp[3] = sf3 * prey*2;
  tmp[4] = sf4 * prez;
  tmp[5] = sf5 * prez*2;

  // communicate values with other procs

  MPI_Allreduce(tmp,sf_coeff,6,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   run the regular toplevel compute method from plain PPPMCG
   which will have individual methods replaced by our threaded
   versions and then call the obligatory force reduction.
------------------------------------------------------------------------- */

void PPPMCGOMP::compute(int eflag, int vflag)
{

  PPPMCG::compute(eflag,vflag);

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMCGOMP::make_rho()
{

  // clear 3d density array

  FFT_SCALAR * _noalias const d = &(density_brick[nzlo_out][nylo_out][nxlo_out]);
  memset(d,0,ngrid*sizeof(FFT_SCALAR));

  // no local atoms with a charge => nothing else to do

  if (num_charged == 0) return;

  const int ix = nxhi_out - nxlo_out + 1;
  const int iy = nyhi_out - nylo_out + 1;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    const double * _noalias const q = atom->q;
    const auto * _noalias const x = (dbl3_t *) atom->x[0];
    const auto * _noalias const p2g = (int3_t *) part2grid[0];

    const double boxlox = boxlo[0];
    const double boxloy = boxlo[1];
    const double boxloz = boxlo[2];

    // determine range of grid points handled by this thread
    int i,jfrom,jto,tid;
    loop_setup_thr(jfrom,jto,tid,ngrid,comm->nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    // loop over my charges, add their contribution to nearby grid points
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // (dx,dy,dz) = distance to "lower left" grid pt

    // loop over all local atoms for all threads

    for (int j = 0; j < num_charged; j++) {
      i = is_charged[j];

      const int nx = p2g[i].a;
      const int ny = p2g[i].b;
      const int nz = p2g[i].t;

      // pre-screen whether this atom will ever come within
      // reach of the data segement this thread is updating.
      if ( ((nz+nlower-nzlo_out)*ix*iy >= jto)
           || ((nz+nupper-nzlo_out+1)*ix*iy < jfrom) ) continue;

      const FFT_SCALAR dx = nx+shiftone - (x[i].x-boxlox)*delxinv;
      const FFT_SCALAR dy = ny+shiftone - (x[i].y-boxloy)*delyinv;
      const FFT_SCALAR dz = nz+shiftone - (x[i].z-boxloz)*delzinv;

      compute_rho1d_thr(r1d,dx,dy,dz);

      const FFT_SCALAR z0 = delvolinv * q[i];

      for (int n = nlower; n <= nupper; ++n) {
        const int jn = (nz+n-nzlo_out)*ix*iy;
        const FFT_SCALAR y0 = z0*r1d[2][n];

        for (int m = nlower; m <= nupper; ++m) {
          const int jm = jn+(ny+m-nylo_out)*ix;
          const FFT_SCALAR x0 = y0*r1d[1][m];

          for (int l = nlower; l <= nupper; ++l) {
            const int jl = jm+nx+l-nxlo_out;
            // make sure each thread only updates
            // "his" elements of the density grid
            if (jl >= jto) break;
            if (jl < jfrom) continue;

            d[jl] += x0*r1d[0][l];
          }
        }
      }
    }
    thr->timer(Timer::KSPACE);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMCGOMP::fieldforce_ik()
{
  // no local atoms with a charge => nothing to do

  if (num_charged == 0) return;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const double * _noalias const q = atom->q;
  const double qqrd2e = force->qqrd2e;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    FFT_SCALAR dx,dy,dz,x0,y0,z0,ekx,eky,ekz;
    int i,ifrom,ito,tid,l,m,n,nx,ny,nz,mx,my,mz;

    loop_setup_thr(ifrom,ito,tid,num_charged,nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    for (int j = ifrom; j < ito; ++j) {
      i = is_charged[j];

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i].x-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i].y-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i].z-boxlo[2])*delzinv;

      compute_rho1d_thr(r1d,dx,dy,dz);

      ekx = eky = ekz = ZEROF;
      for (n = nlower; n <= nupper; n++) {
        mz = n+nz;
        z0 = r1d[2][n];
        for (m = nlower; m <= nupper; m++) {
          my = m+ny;
          y0 = z0*r1d[1][m];
          for (l = nlower; l <= nupper; l++) {
            mx = l+nx;
            x0 = y0*r1d[0][l];
            ekx -= x0*vdx_brick[mz][my][mx];
            eky -= x0*vdy_brick[mz][my][mx];
            ekz -= x0*vdz_brick[mz][my][mx];
          }
        }
      }

      // convert E-field to force

      const double qfactor = qqrd2e * scale * q[i];
      f[i].x += qfactor*ekx;
      f[i].y += qfactor*eky;
      if (slabflag != 2) f[i].z += qfactor*ekz;
    }
    thr->timer(Timer::KSPACE);
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMCGOMP::fieldforce_ad()
{
  // no local atoms with a charge => nothing to do

  if (num_charged == 0) return;

  const double *prd = (triclinic == 0) ? domain->prd : domain->prd_lamda;
  const double hx_inv = nx_pppm/prd[0];
  const double hy_inv = ny_pppm/prd[1];
  const double hz_inv = nz_pppm/prd[2];

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const double * _noalias const q = atom->q;
  const double qqrd2e = force->qqrd2e;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    int i,ifrom,ito,tid,l,m,n,nx,ny,nz,mx,my,mz;
    FFT_SCALAR dx,dy,dz,ekx,eky,ekz;
    double s1,s2,s3,sf;

    loop_setup_thr(ifrom,ito,tid,num_charged,nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());
    FFT_SCALAR * const * const d1d = static_cast<FFT_SCALAR **>(thr->get_drho1d());

    for (int j = ifrom; j < ito; ++j) {
      i = is_charged[j];

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i].x-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i].y-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i].z-boxlo[2])*delzinv;

      compute_rho1d_thr(r1d,dx,dy,dz);
      compute_drho1d_thr(d1d,dx,dy,dz);

      ekx = eky = ekz = ZEROF;
      for (n = nlower; n <= nupper; n++) {
        mz = n+nz;
        for (m = nlower; m <= nupper; m++) {
          my = m+ny;
          for (l = nlower; l <= nupper; l++) {
            mx = l+nx;
            ekx += d1d[0][l]*r1d[1][m]*r1d[2][n]*u_brick[mz][my][mx];
            eky += r1d[0][l]*d1d[1][m]*r1d[2][n]*u_brick[mz][my][mx];
            ekz += r1d[0][l]*r1d[1][m]*d1d[2][n]*u_brick[mz][my][mx];
          }
        }
      }
      ekx *= hx_inv;
      eky *= hy_inv;
      ekz *= hz_inv;

      // convert E-field to force and subtract self forces

      const double qi = q[i];
      const double qfactor = qqrd2e * scale * qi;

      s1 = x[i].x*hx_inv;
      sf = sf_coeff[0]*sin(MY_2PI*s1);
      sf += sf_coeff[1]*sin(MY_4PI*s1);
      sf *= 2.0*qi;
      f[i].x += qfactor*(ekx - sf);

      s2 = x[i].y*hy_inv;
      sf = sf_coeff[2]*sin(MY_2PI*s2);
      sf += sf_coeff[3]*sin(MY_4PI*s2);
      sf *= 2*qi;
      f[i].y += qfactor*(eky - sf);

      s3 = x[i].z*hz_inv;
      sf = sf_coeff[4]*sin(MY_2PI*s3);
      sf += sf_coeff[5]*sin(MY_4PI*s3);
      sf *= 2*qi;
      if (slabflag != 2) f[i].z += qfactor*(ekz - sf);
    }
    thr->timer(Timer::KSPACE);
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void PPPMCGOMP::fieldforce_peratom()
{
  // no local atoms with a charge => nothing to do

  if (num_charged == 0) return;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  const double * _noalias const q = atom->q;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    FFT_SCALAR dx,dy,dz,x0,y0,z0;
    FFT_SCALAR u,v0,v1,v2,v3,v4,v5;
    int i,ifrom,ito,tid,l,m,n,nx,ny,nz,mx,my,mz;

    loop_setup_thr(ifrom,ito,tid,num_charged,nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    for (int j=ifrom; j < ito; ++j) {
      i = is_charged[j];

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i].x-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i].y-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i].z-boxlo[2])*delzinv;

      compute_rho1d_thr(r1d,dx,dy,dz);

      u = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
      for (n = nlower; n <= nupper; n++) {
        mz = n+nz;
        z0 = r1d[2][n];
        for (m = nlower; m <= nupper; m++) {
          my = m+ny;
          y0 = z0*r1d[1][m];
          for (l = nlower; l <= nupper; l++) {
            mx = l+nx;
            x0 = y0*r1d[0][l];
            if (eflag_atom) u += x0*u_brick[mz][my][mx];
            if (vflag_atom) {
              v0 += x0*v0_brick[mz][my][mx];
              v1 += x0*v1_brick[mz][my][mx];
              v2 += x0*v2_brick[mz][my][mx];
              v3 += x0*v3_brick[mz][my][mx];
              v4 += x0*v4_brick[mz][my][mx];
              v5 += x0*v5_brick[mz][my][mx];
            }
          }
        }
      }

      const double qi = q[i];
      if (eflag_atom) eatom[i] += qi*u;
      if (vflag_atom) {
        vatom[i][0] += qi*v0;
        vatom[i][1] += qi*v1;
        vatom[i][2] += qi*v2;
        vatom[i][3] += qi*v3;
        vatom[i][4] += qi*v4;
        vatom[i][5] += qi*v5;
      }
    }
    thr->timer(Timer::KSPACE);
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */
void PPPMCGOMP::compute_rho1d_thr(FFT_SCALAR * const * const r1d, const FFT_SCALAR &dx,
                                const FFT_SCALAR &dy, const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-1; l >= 0; l--) {
      r1 = rho_coeff[l][k] + r1*dx;
      r2 = rho_coeff[l][k] + r2*dy;
      r3 = rho_coeff[l][k] + r3*dz;
    }
    r1d[0][k] = r1;
    r1d[1][k] = r2;
    r1d[2][k] = r3;
  }
}

/* ----------------------------------------------------------------------
   charge assignment into drho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */

void PPPMCGOMP::compute_drho1d_thr(FFT_SCALAR * const * const d1d, const FFT_SCALAR &dx,
                              const FFT_SCALAR &dy, const FFT_SCALAR &dz)
{
  int k,l;
  FFT_SCALAR r1,r2,r3;

  for (k = (1-order)/2; k <= order/2; k++) {
    r1 = r2 = r3 = ZEROF;

    for (l = order-2; l >= 0; l--) {
      r1 = drho_coeff[l][k] + r1*dx;
      r2 = drho_coeff[l][k] + r2*dy;
      r3 = drho_coeff[l][k] + r3*dz;
    }
    d1d[0][k] = r1;
    d1d[1][k] = r2;
    d1d[2][k] = r3;
  }
}
