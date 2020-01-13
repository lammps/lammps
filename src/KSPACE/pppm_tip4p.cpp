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

#include "pppm_tip4p.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "error.h"
#include "math_const.h"

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

PPPMTIP4P::PPPMTIP4P(LAMMPS *lmp) : PPPM(lmp)
{
  triclinic_support = 1;
  tip4pflag = 1;
}

/* ---------------------------------------------------------------------- */

void PPPMTIP4P::init()
{
  // TIP4P PPPM requires newton on, b/c it computes forces on ghost atoms

  if (force->newton == 0)
    error->all(FLERR,"Kspace style pppm/tip4p requires newton on");

  PPPM::init();
}

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

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

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
  if (flag_all) error->all(FLERR,"Out of range atoms - cannot compute PPPM");
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
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  double *xi,xM[3];

  // clear 3d density array

  FFT_SCALAR *vec = &density_brick[nzlo_out][nylo_out][nxlo_out];
  for (i = 0; i < ngrid; i++) vec[i] = ZEROF;

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
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMTIP4P::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz;
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

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const double qfactor = qqrd2e * scale * q[i];
    if (type[i] != typeO) {
      f[i][0] += qfactor*ekx;
      f[i][1] += qfactor*eky;
      if (slabflag != 2) f[i][2] += qfactor*ekz;

    } else {
      fx = qfactor * ekx;
      fy = qfactor * eky;
      fz = qfactor * ekz;
      find_M(i,iH1,iH2,xM);

      f[i][0] += fx*(1 - alpha);
      f[i][1] += fy*(1 - alpha);
      if (slabflag != 2) f[i][2] += fz*(1 - alpha);

      f[iH1][0] += 0.5*alpha*fx;
      f[iH1][1] += 0.5*alpha*fy;
      if (slabflag != 2) f[iH1][2] += 0.5*alpha*fz;

      f[iH2][0] += 0.5*alpha*fx;
      f[iH2][1] += 0.5*alpha*fy;
      if (slabflag != 2) f[iH2][2] += 0.5*alpha*fz;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMTIP4P::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz;
  double *xi;
  int iH1,iH2;
  double xM[3];
  double s1,s2,s3;
  double sf;
  double *prd;
  double fx,fy,fz;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd;

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
    compute_drho1d(dx,dy,dz);

    ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }

    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force and substract self forces

    const double qfactor = qqrd2e * scale;

    s1 = xi[0]*hx_inv;
    s2 = xi[1]*hy_inv;
    s3 = xi[2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2.0*q[i]*q[i];
    fx = qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2.0*q[i]*q[i];
    fy = qfactor*(eky*q[i] - sf);

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2.0*q[i]*q[i];
    fz = qfactor*(ekz*q[i] - sf);

    if (type[i] != typeO) {
      f[i][0] += fx;
      f[i][1] += fy;
      if (slabflag != 2) f[i][2] += fz;

    } else {
      find_M(i,iH1,iH2,xM);

      f[i][0] += fx*(1 - alpha);
      f[i][1] += fy*(1 - alpha);
      if (slabflag != 2) f[i][2] += fz*(1 - alpha);

      f[iH1][0] += 0.5*alpha*fx;
      f[iH1][1] += 0.5*alpha*fy;
      if (slabflag != 2) f[iH1][2] += 0.5*alpha*fz;

      f[iH2][0] += 0.5*alpha*fx;
      f[iH2][1] += 0.5*alpha*fy;
      if (slabflag != 2) f[iH2][2] += 0.5*alpha*fz;
    }
  }
}


/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

void PPPMTIP4P::fieldforce_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  double *xi;
  int iH1,iH2;
  double xM[3];
  FFT_SCALAR u_pa,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;

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

    u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) u_pa += x0*u_brick[mz][my][mx];
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

    if (eflag_atom) {
      if (type[i] != typeO) {
        eatom[i] += q[i]*u_pa;
      } else {
        eatom[i] += q[i]*u_pa*(1-alpha);
        eatom[iH1] += q[i]*u_pa*alpha*0.5;
        eatom[iH2] += q[i]*u_pa*alpha*0.5;
      }
    }
    if (vflag_atom) {
      if (type[i] != typeO) {
        vatom[i][0] += v0*q[i];
        vatom[i][1] += v1*q[i];
        vatom[i][2] += v2*q[i];
        vatom[i][3] += v3*q[i];
        vatom[i][4] += v4*q[i];
        vatom[i][5] += v5*q[i];
      } else {
        vatom[i][0] += v0*(1-alpha)*q[i];
        vatom[i][1] += v1*(1-alpha)*q[i];
        vatom[i][2] += v2*(1-alpha)*q[i];
        vatom[i][3] += v3*(1-alpha)*q[i];
        vatom[i][4] += v4*(1-alpha)*q[i];
        vatom[i][5] += v5*(1-alpha)*q[i];
        vatom[iH1][0] += v0*alpha*0.5*q[i];
        vatom[iH1][1] += v1*alpha*0.5*q[i];
        vatom[iH1][2] += v2*alpha*0.5*q[i];
        vatom[iH1][3] += v3*alpha*0.5*q[i];
        vatom[iH1][4] += v4*alpha*0.5*q[i];
        vatom[iH1][5] += v5*alpha*0.5*q[i];
        vatom[iH2][0] += v0*alpha*0.5*q[i];
        vatom[iH2][1] += v1*alpha*0.5*q[i];
        vatom[iH2][2] += v2*alpha*0.5*q[i];
        vatom[iH2][3] += v3*alpha*0.5*q[i];
        vatom[iH2][4] += v4*alpha*0.5*q[i];
        vatom[iH2][5] += v5*alpha*0.5*q[i];
      }
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
  double **x = atom->x;

  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one(FLERR,"TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one(FLERR,"TIP4P hydrogen has incorrect atom type");

  if (triclinic) {

    // need to use custom code to find the closest image for triclinic,
    // since local atoms are in lambda coordinates, but ghosts are not.

    int *sametag = atom->sametag;
    double xo[3],xh1[3],xh2[3],xm[3];
    const int nlocal = atom->nlocal;

    for (int ii = 0; ii < 3; ++ii) {
      xo[ii] = x[i][ii];
      xh1[ii] = x[iH1][ii];
      xh2[ii] = x[iH2][ii];
    }

    if (i < nlocal) domain->lamda2x(x[i],xo);
    if (iH1 < nlocal) domain->lamda2x(x[iH1],xh1);
    if (iH2 < nlocal) domain->lamda2x(x[iH2],xh2);

    double delx = xo[0] - xh1[0];
    double dely = xo[1] - xh1[1];
    double delz = xo[2] - xh1[2];
    double rsqmin = delx*delx + dely*dely + delz*delz;
    double rsq;
    int closest = iH1;

    // no need to run lamda2x() here -> ghost atoms

    while (sametag[iH1] >= 0) {
      iH1 = sametag[iH1];
      delx = xo[0] - x[iH1][0];
      dely = xo[1] - x[iH1][1];
      delz = xo[2] - x[iH1][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH1;
        xh1[0] = x[iH1][0];
        xh1[1] = x[iH1][1];
        xh1[2] = x[iH1][2];
      }
    }
    iH1 = closest;

    closest = iH2;
    delx = xo[0] - xh2[0];
    dely = xo[1] - xh2[1];
    delz = xo[2] - xh2[2];
    rsqmin = delx*delx + dely*dely + delz*delz;

    while (sametag[iH2] >= 0) {
      iH2 = sametag[iH2];
      delx = xo[0] - x[iH2][0];
      dely = xo[1] - x[iH2][1];
      delz = xo[2] - x[iH2][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < rsqmin) {
        rsqmin = rsq;
        closest = iH2;
        xh2[0] = x[iH2][0];
        xh2[1] = x[iH2][1];
        xh2[2] = x[iH2][2];
      }
    }
    iH2 = closest;

    // finally compute M in real coordinates ...

    double delx1 = xh1[0] - xo[0];
    double dely1 = xh1[1] - xo[1];
    double delz1 = xh1[2] - xo[2];

    double delx2 = xh2[0] - xo[0];
    double dely2 = xh2[1] - xo[1];
    double delz2 = xh2[2] - xo[2];

    xm[0] = xo[0] + alpha * 0.5 * (delx1 + delx2);
    xm[1] = xo[1] + alpha * 0.5 * (dely1 + dely2);
    xm[2] = xo[2] + alpha * 0.5 * (delz1 + delz2);

    // ... and convert M to lamda space for PPPM

    domain->x2lamda(xm,xM);

  } else {

    // set iH1,iH2 to index of closest image to O

    iH1 = domain->closest_image(i,iH1);
    iH2 = domain->closest_image(i,iH2);

    double delx1 = x[iH1][0] - x[i][0];
    double dely1 = x[iH1][1] - x[i][1];
    double delz1 = x[iH1][2] - x[i][2];

    double delx2 = x[iH2][0] - x[i][0];
    double dely2 = x[iH2][1] - x[i][1];
    double delz2 = x[iH2][2] - x[i][2];

    xM[0] = x[i][0] + alpha * 0.5 * (delx1 + delx2);
    xM[1] = x[i][1] + alpha * 0.5 * (dely1 + dely2);
    xM[2] = x[i][2] + alpha * 0.5 * (delz1 + delz2);
  }
}
