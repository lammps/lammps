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
#include "force.h"
#include "memory.h"
#include "math_const.h"

#include <string.h>
#include <math.h>

#include "suffix.h"
using namespace LAMMPS_NS;
using namespace MathConst;

#define OFFSET 16384

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMTIP4POMP::PPPMTIP4POMP(LAMMPS *lmp, int narg, char **arg) :
  PPPMOMP(lmp, narg, arg)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

void PPPMTIP4POMP::init()
{
  // TIP4P PPPM requires newton on, b/c it computes forces on ghost atoms

  if (force->newton == 0)
    error->all(FLERR,"Kspace style pppm/tip4p/omp requires newton on");

  PPPMOMP::init();
}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array 
------------------------------------------------------------------------- */

void PPPMTIP4POMP::particle_map()
{
  const int * const type = atom->type;
  const double * const * const x = atom->x;
  const int nlocal = atom->nlocal;

  int i, flag = 0;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) reduction(+:flag) schedule(static)
#endif
  for (i = 0; i < nlocal; i++) {
    int nx,ny,nz,iH1,iH2;
    double xM[3];

    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);      
    } else {
      xM[0] = x[i][0];
      xM[1] = x[i][1];
      xM[2] = x[i][2];
    }

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((xM[0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((xM[1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((xM[2]-boxlo[2])*delzinv+shift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
	ny+nlower < nylo_out || ny+nupper > nyhi_out ||
	nz+nlower < nzlo_out || nz+nupper > nzhi_out) flag++;
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
  const double * const q = atom->q;
  const double * const * const x = atom->x;
  const int * const type = atom->type;
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

    int l,m,n,nx,ny,nz,mx,my,mz,iH1,iH2;
    FFT_SCALAR dx,dy,dz,x0,y0,z0;
    double xM[3];

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

	if (type[i] == typeO) {
	  find_M(i,iH1,iH2,xM);      
	} else {
	  xM[0] = x[i][0];
	  xM[1] = x[i][1];
	  xM[2] = x[i][2];
	}

	nx = part2grid[i][0];
	ny = part2grid[i][1];
	nz = part2grid[i][2];
	dx = nx+shiftone - (xM[0]-boxlo[0])*delxinv;
	dy = ny+shiftone - (xM[1]-boxlo[1])*delyinv;
	dz = nz+shiftone - (xM[2]-boxlo[2])*delzinv;

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

void PPPMTIP4POMP::fieldforce()
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  const double * const q = atom->q;
  const double * const * const x = atom->x;
  const int * const type = atom->type;
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
    int iH1,iH2;
    double xM[3], fx,fy,fz;
    double ddotf, rOMx, rOMy, rOMz, f1x, f1y, f1z;

    // this if protects against having more threads than local atoms
    if (ifrom < nlocal) { 
      for (int i = ifrom; i < ito; i++) {

	if (type[i] == typeO) {
	  find_M(i,iH1,iH2,xM);      
	} else {
	  xM[0] = x[i][0];
	  xM[1] = x[i][1];
	  xM[2] = x[i][2];
	}

	nx = part2grid[i][0];
	ny = part2grid[i][1];
	nz = part2grid[i][2];
	dx = nx+shiftone - (xM[0]-boxlo[0])*delxinv;
	dy = ny+shiftone - (xM[1]-boxlo[1])*delyinv;
	dz = nz+shiftone - (xM[2]-boxlo[2])*delzinv;

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

	if (type[i] != typeO) {
	  f[i][0] += qfactor*ekx;
	  f[i][1] += qfactor*eky;
	  f[i][2] += qfactor*ekz;

	} else {
	  fx = qfactor * ekx;
	  fy = qfactor * eky;
	  fz = qfactor * ekz;
	  find_M(i,iH1,iH2,xM);

	  rOMx = xM[0] - x[i][0];
	  rOMy = xM[1] - x[i][1];
	  rOMz = xM[2] - x[i][2];
  
	  ddotf = (rOMx * fx + rOMy * fy + rOMz * fz) / (qdist * qdist);

	  f1x = ddotf * rOMx;
	  f1y = ddotf * rOMy;
	  f1z = ddotf * rOMz;

	  f[i][0] += fx - alpha * (fx - f1x);
	  f[i][1] += fy - alpha * (fy - f1y);
	  f[i][2] += fz - alpha * (fz - f1z);

	  f[iH1][0] += 0.5*alpha*(fx - f1x);
	  f[iH1][1] += 0.5*alpha*(fy - f1y);
	  f[iH1][2] += 0.5*alpha*(fz - f1z);

	  f[iH2][0] += 0.5*alpha*(fx - f1x);
	  f[iH2][1] += 0.5*alpha*(fy - f1y);
	  f[iH2][2] += 0.5*alpha*(fz - f1z);
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
  find 2 H atoms bonded to O atom i
  compute position xM of fictitious charge site for O atom
  also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PPPMTIP4POMP::find_M(int i, int &iH1, int &iH2, double *xM)
{
  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one(FLERR,"TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  const double * const * const x = atom->x; 

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];
  domain->minimum_image(delx2,dely2,delz2);

  xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
}
