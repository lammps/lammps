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

#include "pppm_omp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"

#include <string.h>

using namespace LAMMPS_NS;

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMOMP::PPPMOMP(LAMMPS *lmp, int narg, char **arg) :
  PPPM(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE)
{
}


/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMOMP::allocate()
{
  PPPM::allocate();

  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif

    FFT_SCALAR **rho1d_thr;
    memory->create2d_offset(rho1d_thr,3,-order/2,order/2,"pppm:rho1d_thr");
    ThrData *thr = fix->get_thr(tid);
    thr->init_pppm(static_cast<void *>(rho1d_thr));
  }

  const int nzend = (nzhi_out-nzlo_out+1)*nthreads + nzlo_out -1;

  // reallocate density brick, so it fits our needs
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  memory->create3d_offset(density_brick,nzlo_out,nzend,nylo_out,nyhi_out,
			  nxlo_out,nxhi_out,"pppm:density_brick");
}

// NOTE: special version of reduce_data for FFT_SCALAR data type.
// reduce per thread data into the first part of the data
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on
// multi-dimensional arrays like force stored in this order
// x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with writing to the array .
static void data_reduce_fft(FFT_SCALAR *dall, int nall, int nthreads, int ndim, int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    const int nvals = ndim*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nvals) ? nvals : (ifrom + idelta);

    for (int m = ifrom; m < ito; ++m) {
      for (int n = 1; n < nthreads; ++n) {
	dall[m] += dall[n*nvals + m];
	dall[n*nvals + m] = 0.0;
      }
    }
  }
#else
  // NOOP in non-threaded execution.
  return;
#endif
}

/* ----------------------------------------------------------------------
   free memory that depends on # of K-vectors and order 
------------------------------------------------------------------------- */

void PPPMOMP::deallocate()
{
  PPPM::deallocate();
  for (int i=0; i < comm->nthreads; ++i) {
    ThrData * thr = fix->get_thr(i);
    FFT_SCALAR ** rho1d_thr = static_cast<FFT_SCALAR **>(thr->get_rho1d());
    memory->destroy2d_offset(rho1d_thr,-order/2);
  }
}

void PPPMOMP::compute(int eflag, int vflag)
{

  PPPM::compute(eflag,vflag);

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
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid 
------------------------------------------------------------------------- */

void PPPMOMP::make_rho()
{
  const double * const q = atom->q;
  const double * const * const x = atom->x;
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(none)
  {  
    // each thread works on a fixed chunk of atoms.
    const int tid = omp_get_thread_num();
    const int inum = atom->nlocal;
    const int idelta = 1 + inum/nthreads;
    const int ifrom = tid*idelta;
    const int ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
    const int tid = 0;
    const int ifrom = 0;
    const int ito = atom->nlocal;
#endif

    int l,m,n,nx,ny,nz,mx,my,mz;
    FFT_SCALAR dx,dy,dz,x0,y0,z0;

    // set up clear 3d density array
    const int nzoffs = (nzhi_out-nzlo_out+1)*tid;
    FFT_SCALAR * const * const * const db = &(density_brick[nzoffs]);
    memset(&(db[nzlo_out][nylo_out][nxlo_out]),0,ngrid*sizeof(FFT_SCALAR));

    ThrData *thr = fix->get_thr(tid);
    FFT_SCALAR * const * const r1d = static_cast<FFT_SCALAR **>(thr->get_rho1d());

    // loop over my charges, add their contribution to nearby grid points
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // (dx,dy,dz) = distance to "lower left" grid pt
    // (mx,my,mz) = global coords of moving stencil pt

    for (int i = ifrom; i < ito; i++) {

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

      compute_rho1d_thr(r1d,dx,dy,dz);

      z0 = delvolinv * q[i];
      for (n = nlower; n <= nupper; n++) {
	mz = n+nz;
	y0 = z0*r1d[2][n];
	for (m = nlower; m <= nupper; m++) {
	  my = m+ny;
	  x0 = y0*r1d[1][m];
	  for (l = nlower; l <= nupper; l++) {
	    mx = l+nx;
	    db[mz][my][mx] += x0*r1d[0][l];
	  }
	}
      }
    }
#if defined(_OPENMP)
    // reduce 3d density array
    if (nthreads > 1) {
      sync_threads();
      data_reduce_fft(&(density_brick[nzlo_out][nylo_out][nxlo_out]),ngrid,nthreads,1,tid);
    }
  }
#endif
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
------------------------------------------------------------------------- */

void PPPMOMP::fieldforce()
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const double * const q = atom->q;
  const double * const * const x = atom->x;
  const int nthreads = comm->nthreads;
  const double qqrd2e = force->qqrd2e;

#if defined(_OPENMP)
#pragma omp parallel default(none)
  {  
    // each thread works on a fixed chunk of atoms.
    const int tid = omp_get_thread_num();
    const int inum = atom->nlocal;
    const int idelta = 1 + inum/nthreads;
    const int ifrom = tid*idelta;
    const int ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
    const int ifrom = 0;
    const int ito = atom->nlocal;
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    double * const * const f = thr->get_f();
    FFT_SCALAR * const * const r1d =  static_cast<FFT_SCALAR **>(thr->get_rho1d());
    
    int l,m,n,nx,ny,nz,mx,my,mz;
    FFT_SCALAR dx,dy,dz,x0,y0,z0;
    FFT_SCALAR ekx,eky,ekz;

    for (int i = ifrom; i < ito; i++) {

      nx = part2grid[i][0];
      ny = part2grid[i][1];
      nz = part2grid[i][2];
      dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
      dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
      dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

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
      const double qfactor = qqrd2e*scale*q[i];
      f[i][0] += qfactor*ekx;
      f[i][1] += qfactor*eky;
      f[i][2] += qfactor*ekz;
    }
#if defined(_OPENMP)
  }
#endif
}

/* ----------------------------------------------------------------------
   charge assignment into rho1d
   dx,dy,dz = distance of particle from "lower left" grid point 
------------------------------------------------------------------------- */
void PPPMOMP::compute_rho1d_thr(FFT_SCALAR * const * const r1d, const FFT_SCALAR &dx,
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
