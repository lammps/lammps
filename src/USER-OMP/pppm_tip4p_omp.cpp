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

#include "pppm_tip4p_omp.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_omp.h"
#include "force.h"
#include "memory.h"
#include "math_const.h"
#include "math_special.h"

#include <string.h>
#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#else
#define ZEROF 0.0
#endif

#define EPS_HOC 1.0e-7
#define OFFSET 16384

/* ---------------------------------------------------------------------- */

PPPMTIP4POMP::PPPMTIP4POMP(LAMMPS *lmp, int narg, char **arg) :
  PPPMTIP4P(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE)
{
  triclinic_support = 0;
  suffix_flag |= Suffix::OMP;
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMTIP4POMP::allocate()
{
  PPPMTIP4P::allocate();

#if defined(_OPENMP)
#pragma omp parallel default(none)
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
   free memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMTIP4POMP::deallocate()
{
  PPPMTIP4P::deallocate();

#if defined(_OPENMP)
#pragma omp parallel default(none)
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
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
------------------------------------------------------------------------- */

void PPPMTIP4POMP::compute_gf_ik()
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
#pragma omp parallel default(none)
#endif
  {
    double snx,sny,snz;
    double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
    double sum1,dot1,dot2;
    double numerator,denominator;
    double sqk;

    int k,l,m,nx,ny,nz,kper,lper,mper,n,nfrom,nto,tid;

    loop_setup_thr(nfrom, nto, tid, nfft, comm->nthreads);

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
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
------------------------------------------------------------------------- */

void PPPMTIP4POMP::compute_gf_ad()
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
#pragma omp parallel default(none) reduction(+:sf0,sf1,sf2,sf3,sf4,sf5)
#endif
  {
    double snx,sny,snz,sqk;
    double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
    double numerator,denominator;
    int k,l,m,kper,lper,mper,n,nfrom,nto,tid;

    loop_setup_thr(nfrom, nto, tid, nfft, comm->nthreads);

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
  } // end of paralle region
  
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
   run the regular toplevel compute method from plain PPPM
   which will have individual methods replaced by our threaded
   versions and then call the obligatory force reduction.
------------------------------------------------------------------------- */

void PPPMTIP4POMP::compute(int eflag, int vflag)
{

  PPPMTIP4P::compute(eflag,vflag);

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
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void PPPMTIP4POMP::particle_map()
{
  // no local atoms => nothing to do

  if (atom->nlocal == 0) return;

  const int * _noalias const type = atom->type;
  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  int3_t * _noalias const p2g = (int3_t *) part2grid[0];
  const double boxlox = boxlo[0];
  const double boxloy = boxlo[1];
  const double boxloz = boxlo[2];
  const int nlocal = atom->nlocal;

  if (!isfinite(boxlo[0]) || !isfinite(boxlo[1]) || !isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  int i, flag = 0;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:flag) schedule(static)
#endif
  for (i = 0; i < nlocal; i++) {
    dbl3_t xM;
    int iH1,iH2;

    if (type[i] == typeO) {
      find_M_thr(i,iH1,iH2,xM);
    } else {
      xM = x[i];
    }

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    const int nx = static_cast<int> ((xM.x-boxlox)*delxinv+shift) - OFFSET;
    const int ny = static_cast<int> ((xM.y-boxloy)*delyinv+shift) - OFFSET;
    const int nz = static_cast<int> ((xM.z-boxloz)*delzinv+shift) - OFFSET;

    p2g[i].a = nx;
    p2g[i].b = ny;
    p2g[i].t = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag++;
  }

  int flag_all;
  MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_SUM,world);
  if (flag_all) error->all(FLERR,"Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMTIP4POMP::make_rho()
{

  // clear 3d density array

  FFT_SCALAR * _noalias const d = &(density_brick[nzlo_out][nylo_out][nxlo_out]);
  memset(d,0,ngrid*sizeof(FFT_SCALAR));

  // no local atoms => nothing else to do

  const int nlocal = atom->nlocal;
  if (nlocal == 0) return;

  const int ix = nxhi_out - nxlo_out + 1;
  const int iy = nyhi_out - nylo_out + 1;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    const double * _noalias const q = atom->q;
    const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
    const int3_t * _noalias const p2g = (int3_t *) part2grid[0];
    const int * _noalias const type = atom->type;
    dbl3_t xM;

    const double boxlox = boxlo[0];
    const double boxloy = boxlo[1];
    const double boxloz = boxlo[2];

    // determine range of grid points handled by this thread
    int i,jfrom,jto,tid,iH1,iH2;
    loop_setup_thr(jfrom,jto,tid,ngrid,comm->nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    // loop over my charges, add their contribution to nearby grid points
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // (dx,dy,dz) = distance to "lower left" grid pt

    // loop over all local atoms for all threads
    for (i = 0; i < nlocal; i++) {

      const int nx = p2g[i].a;
      const int ny = p2g[i].b;
      const int nz = p2g[i].t;

      // pre-screen whether this atom will ever come within 
      // reach of the data segement this thread is updating.
      if ( ((nz+nlower-nzlo_out)*ix*iy >= jto)
           || ((nz+nupper-nzlo_out+1)*ix*iy < jfrom) ) continue;

      if (type[i] == typeO) {
        find_M_thr(i,iH1,iH2,xM);
      } else {
        xM = x[i];
      }
      const FFT_SCALAR dx = nx+shiftone - (xM.x-boxlox)*delxinv;
      const FFT_SCALAR dy = ny+shiftone - (xM.y-boxloy)*delyinv;
      const FFT_SCALAR dz = nz+shiftone - (xM.z-boxloz)*delzinv;

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
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMTIP4POMP::fieldforce_ik()
{
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;

  // no local atoms => nothing to do

  if (nlocal == 0) return;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  const double * _noalias const q = atom->q;
  const int3_t * _noalias const p2g = (int3_t *) part2grid[0];
  const int * _noalias const type = atom->type;

  const double qqrd2e = force->qqrd2e;
  const double boxlox = boxlo[0];
  const double boxloy = boxlo[1];
  const double boxloz = boxlo[2];

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    dbl3_t xM;
    FFT_SCALAR x0,y0,z0,ekx,eky,ekz;
    int i,ifrom,ito,tid,iH1,iH2,l,m,n,mx,my,mz;

    loop_setup_thr(ifrom,ito,tid,nlocal,nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    for (i = ifrom; i < ito; ++i) {
      if (type[i] == typeO) {
        find_M_thr(i,iH1,iH2,xM);
      } else xM = x[i];

      const int nx = p2g[i].a;
      const int ny = p2g[i].b;
      const int nz = p2g[i].t;
      const FFT_SCALAR dx = nx+shiftone - (xM.x-boxlox)*delxinv;
      const FFT_SCALAR dy = ny+shiftone - (xM.y-boxloy)*delyinv;
      const FFT_SCALAR dz = nz+shiftone - (xM.z-boxloz)*delzinv;

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
      if (type[i] != typeO) {
        f[i].x += qfactor*ekx;
        f[i].y += qfactor*eky;
        if (slabflag != 2) f[i].z += qfactor*ekz;

      } else {
        const double fx = qfactor * ekx;
        const double fy = qfactor * eky;
        const double fz = qfactor * ekz;

        f[i].x += fx*(1 - alpha);
        f[i].y += fy*(1 - alpha);
        if (slabflag != 2) f[i].z += fz*(1 - alpha);

        f[iH1].x += 0.5*alpha*fx;
        f[iH1].y += 0.5*alpha*fy;
        if (slabflag != 2) f[iH1].z += 0.5*alpha*fz;

        f[iH2].x += 0.5*alpha*fx;
        f[iH2].y += 0.5*alpha*fy;
        if (slabflag != 2) f[iH2].z += 0.5*alpha*fz;
      }
    }
  } // end of parallel region
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMTIP4POMP::fieldforce_ad()
{
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;

  // no local atoms => nothing to do

  if (nlocal == 0) return;

  const double *prd = (triclinic == 0) ? domain->prd : domain->prd_lamda;
  const double hx_inv = nx_pppm/prd[0];
  const double hy_inv = ny_pppm/prd[1];
  const double hz_inv = nz_pppm/prd[2];

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  const double * _noalias const q = atom->q;
  const int3_t * _noalias const p2g = (int3_t *) part2grid[0];
  const int * _noalias const type = atom->type;

  const double qqrd2e = force->qqrd2e;
  const double boxlox = boxlo[0];
  const double boxloy = boxlo[1];
  const double boxloz = boxlo[2];

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    double s1,s2,s3,sf;
    dbl3_t xM;
    FFT_SCALAR ekx,eky,ekz;
    int i,ifrom,ito,tid,iH1,iH2,l,m,n,mx,my,mz;

    loop_setup_thr(ifrom,ito,tid,nlocal,nthreads);

    // get per thread data
    ThrData *thr = fix->get_thr(tid);
    dbl3_t * _noalias const f = (dbl3_t *) thr->get_f()[0];
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());
    FFT_SCALAR * const * const d1d = static_cast<FFT_SCALAR **>(thr->get_drho1d());

    for (i = ifrom; i < ito; ++i) {
      if (type[i] == typeO) {
        find_M_thr(i,iH1,iH2,xM);
      } else xM = x[i];

      const int nx = p2g[i].a;
      const int ny = p2g[i].b;
      const int nz = p2g[i].t;
      const FFT_SCALAR dx = nx+shiftone - (xM.x-boxlox)*delxinv;
      const FFT_SCALAR dy = ny+shiftone - (xM.y-boxloy)*delyinv;
      const FFT_SCALAR dz = nz+shiftone - (xM.z-boxloz)*delzinv;

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

      // convert E-field to force and substract self forces

      const double qi = q[i];
      const double qfactor = qqrd2e * scale * qi;

      s1 = x[i].x*hx_inv;
      sf = sf_coeff[0]*sin(MY_2PI*s1);
      sf += sf_coeff[1]*sin(MY_4PI*s1);
      sf *= 2.0*qi;
      const double fx = qfactor*(ekx - sf);

      s2 = x[i].y*hy_inv;
      sf = sf_coeff[2]*sin(MY_2PI*s2);
      sf += sf_coeff[3]*sin(MY_4PI*s2);
      sf *= 2.0*qi;
      const double fy = qfactor*(eky - sf);

      s3 = x[i].z*hz_inv;
      sf = sf_coeff[4]*sin(MY_2PI*s3);
      sf += sf_coeff[5]*sin(MY_4PI*s3);
      sf *= 2.0*qi;
      const double fz = qfactor*(ekz - sf);

      if (type[i] != typeO) {
        f[i].x += fx;
        f[i].y += fy;
        if (slabflag != 2) f[i].z += fz;

      } else {
        f[i].x += fx*(1 - alpha);
        f[i].y += fy*(1 - alpha);
        if (slabflag != 2) f[i].z += fz*(1 - alpha);

        f[iH1].x += 0.5*alpha*fx;
        f[iH1].y += 0.5*alpha*fy;
        if (slabflag != 2) f[iH1].z += 0.5*alpha*fz;

        f[iH2].x += 0.5*alpha*fx;
        f[iH2].y += 0.5*alpha*fy;
        if (slabflag != 2) f[iH2].z += 0.5*alpha*fz;
      }
    }
  } // end of parallel region
}

/* ----------------------------------------------------------------------
  find 2 H atoms bonded to O atom i
  compute position xM of fictitious charge site for O atom
  also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PPPMTIP4POMP::find_M_thr(int i, int &iH1, int &iH2, dbl3_t &xM)
{
  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one(FLERR,"TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];

  double delx1 = x[iH1].x - x[i].x;
  double dely1 = x[iH1].y - x[i].y;
  double delz1 = x[iH1].z - x[i].z;
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = x[iH2].x - x[i].x;
  double dely2 = x[iH2].y - x[i].y;
  double delz2 = x[iH2].z - x[i].z;
  domain->minimum_image(delx2,dely2,delz2);

  xM.x = x[i].x + alpha * 0.5 * (delx1 + delx2);
  xM.y = x[i].y + alpha * 0.5 * (dely1 + dely2);
  xM.z = x[i].z + alpha * 0.5 * (delz1 + delz2);
}


/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point
------------------------------------------------------------------------- */
void PPPMTIP4POMP::compute_rho1d_thr(FFT_SCALAR * const * const r1d, const FFT_SCALAR &dx,
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

void PPPMTIP4POMP::compute_drho1d_thr(FFT_SCALAR * const * const d1d, const FFT_SCALAR &dx,
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
