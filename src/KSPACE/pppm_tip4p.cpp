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
   Contributing authors: Amalie Frischknecht and Ahmed Ismail (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "pppm_tip4p.h"
#include "atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define OFFSET 4096

/* ---------------------------------------------------------------------- */

PPPMTIP4P::PPPMTIP4P(LAMMPS *lmp, int narg, char **arg) :
  PPPM(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array 
------------------------------------------------------------------------- */

void PPPMTIP4P::particle_map()
{
  int nx,ny,nz,iH1,iH2;
  double *xi,xM[3];

  int *type = atom->type;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);      
      xi = xM;
    } else xi = x[i];

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((xi[0]-boxlo[0])*delxinv+shift) - OFFSET;
    ny = static_cast<int> ((xi[1]-boxlo[1])*delyinv+shift) - OFFSET;
    nz = static_cast<int> ((xi[2]-boxlo[2])*delzinv+shift) - OFFSET;

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
  if (flag_all) error->all("Out of range atoms - cannot compute PPPM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid 
------------------------------------------------------------------------- */

void PPPMTIP4P::make_rho()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz,iH1,iH2;
  double dx,dy,dz,x0,y0,z0;
  double *xi,xM[3];

  // clear 3d density array

  double *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = 0.0;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int *type = atom->type; 
  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);      
      xi = xM;
    } else xi = x[i];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (xi[0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (xi[1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (xi[2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	x0 = y0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  density_brick[mz][my][mx] += x0*rho1d[0][l];
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles 
------------------------------------------------------------------------- */

void PPPMTIP4P::fieldforce()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;
  double ek[3];
  double *xi;
  int iH1,iH2;
  double xM[3];
  double fx,fy,fz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    if (type[i] == typeO) {
      find_M(i,iH1,iH2,xM);      
      xi = xM;
    } else xi = x[i];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (xi[0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (xi[1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (xi[2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    ek[0] = ek[1] = ek[2] = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
	my = m+ny;
	y0 = z0*rho1d[1][m];
	for (l = nlower; l <= nupper; l++) {
	  mx = l+nx;
	  x0 = y0*rho1d[0][l];
	  ek[0] -= x0*vdx_brick[mz][my][mx];
	  ek[1] -= x0*vdy_brick[mz][my][mx];
	  ek[2] -= x0*vdz_brick[mz][my][mx];
	}
      }
    }

    // convert E-field to force

    if (type[i] != typeO) {
      f[i][0] += qqrd2e*q[i]*ek[0];
      f[i][1] += qqrd2e*q[i]*ek[1];
      f[i][2] += qqrd2e*q[i]*ek[2];
    } else {

      fx = qqrd2e * q[i] * ek[0];
      fy = qqrd2e * q[i] * ek[1];
      fz = qqrd2e * q[i] * ek[2];
      find_M(i,iH1,iH2,xM);

      f[i][0] += fx*(1.0-2.0*alpha);
      f[i][1] += fy*(1.0-2.0*alpha);
      f[i][2] += fz*(1.0-2.0*alpha);

      f[iH1][0] += alpha*(fx); 
      f[iH1][1] += alpha*(fy); 
      f[iH1][2] += alpha*(fz); 

      f[iH2][0] += alpha*(fx); 
      f[iH2][1] += alpha*(fy); 
      f[iH2][2] += alpha*(fz); 
    }
  }
}

/* ----------------------------------------------------------------------
  find 2 H atoms bonded to O atom i
  compute position xM of fictitious charge site for O atom
  also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PPPMTIP4P::find_M(int i, int &iH1, int &iH2, double *xM)
{
  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one("TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one("TIP4P hydrogen has incorrect atom type");

  double **x = atom->x; 

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];
  domain->minimum_image(delx2,dely2,delz2);

  xM[0] = x[i][0] + alpha * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * (delz1 + delz2);
}
