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
#include "domain.h"
#include "force.h"
#include "memory.h"
#include "math_const.h"

#include <string.h>
#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

#define EPS_HOC 1.0e-7

/* ---------------------------------------------------------------------- */

PPPMOMP::PPPMOMP(LAMMPS *lmp, int narg, char **arg) :
  PPPM(lmp, narg, arg), ThrOMP(lmp, THR_KSPACE)
{
  suffix_flag |= Suffix::OMP;
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
// writing to the array.
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

    // this if protects against having more threads than atoms
    if (ifrom < nall) { 
      for (int m = ifrom; m < ito; ++m) {
	for (int n = 1; n < nthreads; ++n) {
	  dall[m] += dall[n*nvals + m];
	  dall[n*nvals + m] = 0.0;
	}
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

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void PPPMOMP::setup()
{
  int i,j,k,n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;
    
  delxinv = nx_pppm/xprd;
  delyinv = ny_pppm/yprd;
  delzinv = nz_pppm/zprd_slab;

  delvolinv = delxinv*delyinv*delzinv;

  const double unitkx = (2.0*MY_PI/xprd);
  const double unitky = (2.0*MY_PI/yprd);
  const double unitkz = (2.0*MY_PI/zprd_slab);

  // fkx,fky,fkz for my FFT grid pts

  double per;

  for (i = nxlo_fft; i <= nxhi_fft; i++) {
    per = i - nx_pppm*(2*i/nx_pppm);
    fkx[i] = unitkx*per;
  }

  for (i = nylo_fft; i <= nyhi_fft; i++) {
    per = i - ny_pppm*(2*i/ny_pppm);
    fky[i] = unitky*per;
  }

  for (i = nzlo_fft; i <= nzhi_fft; i++) {
    per = i - nz_pppm*(2*i/nz_pppm);
    fkz[i] = unitkz*per;
  }

  // virial coefficients

  double sqk,vterm;

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
	sqk = fkx[i]*fkx[i] + fky[j]*fky[j] + fkz[k]*fkz[k];
	if (sqk == 0.0) {
	  vg[n][0] = 0.0;
	  vg[n][1] = 0.0;
	  vg[n][2] = 0.0;
	  vg[n][3] = 0.0;
	  vg[n][4] = 0.0;
	  vg[n][5] = 0.0;
	} else {
	  vterm = -2.0 * (1.0/sqk + 0.25/(g_ewald*g_ewald));
	  vg[n][0] = 1.0 + vterm*fkx[i]*fkx[i];
	  vg[n][1] = 1.0 + vterm*fky[j]*fky[j];
	  vg[n][2] = 1.0 + vterm*fkz[k]*fkz[k];
	  vg[n][3] = vterm*fkx[i]*fky[j];
	  vg[n][4] = vterm*fkx[i]*fkz[k];
	  vg[n][5] = vterm*fky[j]*fkz[k];
	}
	n++;
      }
    }
  }

  // modified (Hockney-Eastwood) Coulomb Green's function

  const int nbx = static_cast<int> ((g_ewald*xprd/(MY_PI*nx_pppm)) * 
				    pow(-log(EPS_HOC),0.25));
  const int nby = static_cast<int> ((g_ewald*yprd/(MY_PI*ny_pppm)) * 
				    pow(-log(EPS_HOC),0.25));
  const int nbz = static_cast<int> ((g_ewald*zprd_slab/(MY_PI*nz_pppm)) * 
				    pow(-log(EPS_HOC),0.25));

  const double form = 1.0;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
    int tid,nn,nnfrom,nnto,nx,ny,nz,k,l,m;
    double snx,sny,snz,sqk;
    double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
    double sum1,dot1,dot2;
    double numerator,denominator;
    const double gew2 = -4.0*g_ewald*g_ewald;
    
    const int nnx = nxhi_fft-nxlo_fft+1;
    const int nny = nyhi_fft-nylo_fft+1;
    
    loop_setup_thr(nnfrom, nnto, tid, nfft, comm->nthreads);
    
    for (m = nzlo_fft; m <= nzhi_fft; m++) {
      
      const double fkzm = fkz[m];
      snz = sin(0.5*fkzm*zprd_slab/nz_pppm);
      snz *= snz;

      for (l = nylo_fft; l <= nyhi_fft; l++) {
	const double fkyl = fky[l];
	sny = sin(0.5*fkyl*yprd/ny_pppm);
	sny *= sny;

	for (k = nxlo_fft; k <= nxhi_fft; k++) {

	  /* only compute the part designated to this thread */
	  nn = k-nxlo_fft + nnx*(l-nylo_fft + nny*(m-nzlo_fft));
	  if ((nn < nnfrom) || (nn >=nnto)) continue;

	  const double fkxk = fkx[k];
	  snx = sin(0.5*fkxk*xprd/nx_pppm);
	  snx *= snx;
      
	  sqk = fkxk*fkxk + fkyl*fkyl + fkzm*fkzm;

	  if (sqk != 0.0) {
	    numerator = form*MY_4PI/sqk;
	    denominator = gf_denom(snx,sny,snz);  
	    sum1 = 0.0;
	    for (nx = -nbx; nx <= nbx; nx++) {
	      qx = fkxk + unitkx*nx_pppm*nx;
	      sx = exp(qx*qx/gew2);
	      wx = 1.0;
	      argx = 0.5*qx*xprd/nx_pppm;
	      if (argx != 0.0) wx = pow(sin(argx)/argx,order);
	      wx *=wx;

	      for (ny = -nby; ny <= nby; ny++) {
		qy = fkyl + unitky*ny_pppm*ny;
		sy = exp(qy*qy/gew2);
		wy = 1.0;
		argy = 0.5*qy*yprd/ny_pppm;
		if (argy != 0.0) wy = pow(sin(argy)/argy,order);
		wy *= wy;

		for (nz = -nbz; nz <= nbz; nz++) {
		  qz = fkzm + unitkz*nz_pppm*nz;
		  sz = exp(qz*qz/gew2);
		  wz = 1.0;
		  argz = 0.5*qz*zprd_slab/nz_pppm;
		  if (argz != 0.0) wz = pow(sin(argz)/argz,order);
		  wz *= wz;

		  dot1 = fkxk*qx + fkyl*qy + fkzm*qz;
		  dot2 = qx*qx+qy*qy+qz*qz;
		  sum1 += (dot1/dot2) * sx*sy*sz * wx*wy*wz;
		}
	      }
	    }
	    greensfn[nn] = numerator*sum1/denominator;
	  } else greensfn[nn] = 0.0;
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   run the regular toplevel compute method from plain PPPPM 
   which will have individual methods replaced by our threaded
   versions and then call the obligatory force reduction.
------------------------------------------------------------------------- */

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
  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {  
#if defined(_OPENMP)
    // each thread works on a fixed chunk of atoms.
    const int tid = omp_get_thread_num();
    const int inum = nlocal;
    const int idelta = 1 + inum/nthreads;
    const int ifrom = tid*idelta;
    const int ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
    const int tid = 0;
    const int ifrom = 0;
    const int ito = nlocal;
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
    
    // this if protects against having more threads than local atoms
    if (ifrom < nlocal) { 
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
    }
#if defined(_OPENMP)
    // reduce 3d density array
    if (nthreads > 1) {
      data_reduce_fft(&(density_brick[nzlo_out][nylo_out][nxlo_out]),ngrid,nthreads,1,tid);
    }
#endif
  }
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
  const int nlocal = atom->nlocal;
  const double qqrd2e = force->qqrd2e;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {  
#if defined(_OPENMP)
    // each thread works on a fixed chunk of atoms.
    const int tid = omp_get_thread_num();
    const int inum = nlocal;
    const int idelta = 1 + inum/nthreads;
    const int ifrom = tid*idelta;
    const int ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
    const int ifrom = 0;
    const int ito = nlocal;
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    double * const * const f = thr->get_f();
    FFT_SCALAR * const * const r1d =  static_cast<FFT_SCALAR **>(thr->get_rho1d());
    
    int l,m,n,nx,ny,nz,mx,my,mz;
    FFT_SCALAR dx,dy,dz,x0,y0,z0;
    FFT_SCALAR ekx,eky,ekz;

    // this if protects against having more threads than local atoms
    if (ifrom < nlocal) { 
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
    }
  }
}

/* ----------------------------------------------------------------------
 interpolate from grid to get per-atom energy/virial
 ------------------------------------------------------------------------- */

void PPPMOMP::fieldforce_peratom()
{
  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  const double * const q = atom->q;
  const double * const * const x = atom->x;
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;

#if defined(_OPENMP)
#pragma omp parallel default(none)
#endif
  {
#if defined(_OPENMP)
    // each thread works on a fixed chunk of atoms.
    const int tid = omp_get_thread_num();
    const int inum = nlocal;
    const int idelta = 1 + inum/nthreads;
    const int ifrom = tid*idelta;
    const int ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;
#else
    const int ifrom = 0;
    const int ito = nlocal;
    const int tid = 0;
#endif
    ThrData *thr = fix->get_thr(tid);
    FFT_SCALAR * const * const r1d =  static_cast<FFT_SCALAR **>(thr->get_rho1d());

    int i,l,m,n,nx,ny,nz,mx,my,mz;
    FFT_SCALAR dx,dy,dz,x0,y0,z0;
    FFT_SCALAR u,v0,v1,v2,v3,v4,v5;

    // this if protects against having more threads than local atoms
    if (ifrom < nlocal) {
      for (int i = ifrom; i < ito; i++) {

	nx = part2grid[i][0];
	ny = part2grid[i][1];
	nz = part2grid[i][2];
	dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
	dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
	dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

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

	if (eflag_atom) eatom[i] += q[i]*u;
	if (vflag_atom) {
	  vatom[i][0] += v0;
	  vatom[i][1] += v1;
	  vatom[i][2] += v2;
	  vatom[i][3] += v3;
	  vatom[i][4] += v4;
	  vatom[i][5] += v5;
	}
      }
    }
  }
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
