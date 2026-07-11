// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"

#include "atom.h"
#include "domain.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

using MathSpecial::powint;

#define ANINT(x) ((x)>0 ? floor((x)+0.5) : ceil((x)-0.5))

/* ----------------------------------------------------------------------
   lattice = setup periodic boundary conditions
   lattice stores the periodic box dimensions and sets angle
   values to be used in computing fractional coordinates
------------------------------------------------------------------------- */

void PairAmoeba::lattice()
{
  recip[0][0] = recip[0][1] = recip[0][2] = 0.0;
  recip[1][0] = recip[1][1] = recip[1][2] = 0.0;
  recip[2][0] = recip[2][1] = recip[2][2] = 0.0;

  double *h_inv = domain->h_inv;

  recip[0][0] = h_inv[0];
  recip[1][1] = h_inv[1];
  recip[2][2] = h_inv[2];

  if (domain->triclinic) {
    recip[1][0] = h_inv[5];
    recip[2][0] = h_inv[4];
    recip[2][1] = h_inv[3];
  }
}

/* ----------------------------------------------------------------------
   moduli = store the inverse DFT moduli
   moduli sets the moduli of the inverse discrete Fourier transform of the B-splines
------------------------------------------------------------------------- */

void PairAmoeba::moduli()
{
  int i;

  // perform dynamic allocation of local arrays

  int maxfft = MAX(nfft1,nfft2);
  maxfft = MAX(maxfft,nfft3);

  if (maxfft > _nfft_max) {
    memory->destroy(_moduli_bsarray);
    _nfft_max = maxfft;
    memory->create(_moduli_bsarray,_nfft_max,"amoeba:_moduli_bsarray");
  }

  // compute and load the moduli values

  double x = 0.0;
  bspline(x,bsorder,_moduli_array);

  for (i = 0; i < maxfft; i++) _moduli_bsarray[i] = 0.0;
  for (i = 0; i < bsorder; i++) _moduli_bsarray[i+1] = _moduli_array[i];

  dftmod(bsmod1,_moduli_bsarray,nfft1,bsorder);
  dftmod(bsmod2,_moduli_bsarray,nfft2,bsorder);
  dftmod(bsmod3,_moduli_bsarray,nfft3,bsorder);
}

/* ----------------------------------------------------------------------
   bspline = determine B-spline coefficients
   bspline calculates the coefficients for an n-th order B-spline approximation
------------------------------------------------------------------------- */

void PairAmoeba::bspline(double x, int n, double *c)
{
  int i,k;
  double denom;

  // initialize the B-spline as the linear case

  c[0] = 1.0 - x;
  c[1] = x;

  // compute standard B-spline recursion to n-th order

  /*for (k = 2; k < n; k++) {
    denom = 1.0 / (k-1);
    c[k] = x * c[k-1] * denom;
    for (i = 0; i < k-1; i++)
      c[k-i] = ((x+i)*c[k-i-1] + (k-i-x)*c[k-i]) * denom;
    c[0] = (1.0-x) * c[0] * denom;
  }*/

  for (k = 3; k <= n; k++) {
    denom = 1.0 / (k-1);
    c[k-1] = x * c[k-2] * denom;
    for (i = 1; i <= k-2; i++)
      c[k-i-1] = ((x+i)*c[k-i-2] + (k-i-x)*c[k-i-1]) * denom;
    c[0] = (1.0-x) * c[0] * denom;
  }

}

/* ----------------------------------------------------------------------
   dftmod = discrete Fourier transform modulus
   dftmod computes the modulus of the discrete Fourier transform
     of bsarray and stores it in bsmod
------------------------------------------------------------------------- */

void PairAmoeba::dftmod(double *bsmod, double *bsarray, int nfft, int order)
{
  int i,j,k;
  int jcut;
  int order2;
  double eps,zeta;
  double arg,factor;
  double sum1,sum2;

  // get the modulus of the discrete Fourier transform

  factor = 2.0 * MY_PI / nfft;
  for (i = 0; i < nfft; i++) {
    sum1 = 0.0;
    sum2 = 0.0;
    for (j = 0; j < nfft; j++) {
      arg = factor*i*j;
      sum1 += bsarray[j]*cos(arg);
      sum2 += bsarray[j]*sin(arg);
    }
    bsmod[i] = sum1*sum1 + sum2*sum2;
    //printf("BSMOD AAA i %d bsmod %g\n",i,bsmod[i]);
  }

  // fix for exponential Euler spline interpolation failure

  eps = 1.0e-7;
  if (bsmod[0] < eps) bsmod[0] = 0.5 * bsmod[1];
  for (i = 1; i < nfft-1; i++)
    if (bsmod[i] < eps) bsmod[i] = 0.5 * (bsmod[i-1] + bsmod[i+1]);
  if (bsmod[nfft-1] < eps) bsmod[nfft-1] = 0.5 * bsmod[nfft-2];

  // compute and apply the optimal zeta coefficient

  jcut = 50;
  order2 = 2 * order;
  for (i = 1; i <= nfft; i++) {
    k = i - 1;
    if (i > nfft/2) k -= nfft;
    if (k == 0) zeta = 1.0;
    else {
      sum1 = 1.0;
      sum2 = 1.0;
      factor = MY_PI * k / nfft;
      for (j = 1; j <= jcut; j++) {
        arg = factor / (factor + MY_PI*j);
        sum1 += powint(arg,order);
        sum2 += powint(arg,order2);
      }
      for (j = 1; j <= jcut; j++) {
        arg = factor / (factor - MY_PI*j);
        sum1 += powint(arg,order);
        sum2 += powint(arg,order2);
      }
      zeta = sum2 / sum1;
    }
    bsmod[i-1] *= zeta*zeta;
  }
}

/* ----------------------------------------------------------------------
   bspline_fill = get PME B-spline coefficients
   bspline_fill finds B-spline coefficients and derivatives
   for PME atomic sites along the fractional coordinate axes
------------------------------------------------------------------------- */

void PairAmoeba::bspline_fill()
{
  int ifr;
  double w,fr,eps;
  double lamda[3];

  int nlocal = atom->nlocal;
  double **x = atom->x;

  // offset used to shift sites off exact lattice bounds

  eps = 1.0e-8;

  // get the B-spline coefficients for each atomic site

  for (int i = 0; i < nlocal; i++) {

    // NOTE: what about offset/shift and w < 0 or w > 1
    // NOTE: could subtract off nlpts to start with
    // NOTE: this is place to check that stencil size does not
    //       go out of bounds relative to igrid for a proc's sub-domain
    //       similar to PPPM::particle_map()
    //       subtracting eps is strange, and could mess up the check
    //       better to check here than in methods like grid_mpole()
    //       but check needs to be valid for all KSpace terms
    // NOTE: could convert x -> lamda for entire set of Nlocal atoms

    domain->x2lamda(x[i],lamda);

    w = lamda[0];
    fr = nfft1 * w;
    ifr = static_cast<int> (fr-eps);
    w = fr - ifr;
    igrid[i][0] = ifr;
    //igrid[i][0] = ifr + 1;
    //if (igrid[i][0] == nfft1) igrid[i][0] = 0;
    bsplgen(w,thetai1[i]);

    w = lamda[1];
    fr = nfft2 * w;
    ifr = static_cast<int> (fr-eps);
    w = fr - ifr;
    igrid[i][1] = ifr;
    //igrid[i][1] = ifr + 1;
    //if (igrid[i][1] == nfft2) igrid[i][1] = 0;
    bsplgen(w,thetai2[i]);

    w = lamda[2];
    fr = nfft3 * w;
    ifr = static_cast<int> (fr-eps);
    w = fr - ifr;
    igrid[i][2] = ifr;
    //igrid[i][2] = ifr + 1;
    //if (igrid[i][2] == nfft3) igrid[i][2] = 0;
    bsplgen(w,thetai3[i]);
  }
}

/* ----------------------------------------------------------------------
   bsplgen = B-spline coefficients for an atom
   bsplgen gets B-spline coefficients and derivatives for
   a single PME atomic site along a particular direction
------------------------------------------------------------------------- */

void PairAmoeba::bsplgen(double w, double **thetai)
{
  int i,j,k;
  int level;
  double denom;

  level = 4;

  // initialization to get to 2nd order recursion

  bsbuild[1][1] = w;
  bsbuild[0][1] = 1.0 - w;

  // perform one pass to get to 3rd order recursion

  bsbuild[2][2] = 0.5 * w * bsbuild[1][1];
  bsbuild[1][2] = 0.5 * ((1.0+w)*bsbuild[0][1] + (2.0-w)*bsbuild[1][1]);
  bsbuild[0][2] = 0.5 * (1.0-w) * bsbuild[0][1];

  // compute standard B-spline recursion to desired order

  for (i = 4; i <= bsorder; i++) {
    k = i - 1;
    denom = 1.0 / k;
    bsbuild[i-1][i-1] = denom * w * bsbuild[k-1][k-1];
    for (j = 1; j <= i-2; j++) {
      bsbuild[i-j-1][i-1] = denom *
        ((w+j)*bsbuild[i-j-2][k-1] + (i-j-w)*bsbuild[i-j-1][k-1]);
    }
    bsbuild[0][i-1] = denom * (1.0-w) * bsbuild[0][k-1];
  }

  // get coefficients for the B-spline first derivative

  k = bsorder - 1;
  bsbuild[bsorder-1][k-1] = bsbuild[bsorder-2][k-1];
  for (int i = bsorder-1; i >= 2; i--)
    bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
  bsbuild[0][k-1] = -bsbuild[0][k-1];

  // get coefficients for the B-spline second derivative

  if (level == 4) {
    k = bsorder - 2;
    bsbuild[bsorder-2][k-1] = bsbuild[bsorder-3][k-1];
    for (int i = bsorder-2; i >= 2; i--)
      bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
    bsbuild[0][k-1] = -bsbuild[0][k-1];
    bsbuild[bsorder-1][k-1] = bsbuild[bsorder-2][k-1];
    for (int i = bsorder-1; i >= 2; i--)
      bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
    bsbuild[0][k-1] = -bsbuild[0][k-1];

    // get coefficients for the B-spline third derivative

    k = bsorder - 3;
    bsbuild[bsorder-3][k-1] = bsbuild[bsorder-4][k-1];
    for (int i = bsorder-3; i >= 2; i--)
      bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
    bsbuild[0][k-1] = -bsbuild[0][k-1];
    bsbuild[bsorder-2][k-1] = bsbuild[bsorder-3][k-1];
    for (int i = bsorder-2; i >= 2; i--)
      bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
    bsbuild[0][k-1] = -bsbuild[0][k-1];
    bsbuild[bsorder-1][k-1] = bsbuild[bsorder-2][k-1];
    for (int i = bsorder-1; i >= 2; i--)
      bsbuild[i-1][k-1] = bsbuild[i-2][k-1] - bsbuild[i-1][k-1];
    bsbuild[0][k-1] = -bsbuild[0][k-1];
  }

  // copy coefficients from temporary to permanent storage

  for (int i = 1; i <= bsorder; i++)
    for (int j = 1; j <= level; j++)
      thetai[i-1][j-1] = bsbuild[i-1][bsorder-j];
}

/* ----------------------------------------------------------------------
   cmp_to_fmp = transformation of multipoles
   cmp_to_fmp transforms the atomic multipoles from Cartesian
   to fractional coordinates
------------------------------------------------------------------------- */

void PairAmoeba::cmp_to_fmp(double **cmp, double **fmp)
{
  int i,j,k;

  // find the matrix to convert Cartesian to fractional

  cart_to_frac();

  // apply the transformation to get the fractional multipoles

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    fmp[i][0] = ctf[0][0] * cmp[i][0];
    for (j = 1; j < 4; j++) {
      fmp[i][j] = 0.0;
      for (k = 1; k < 4; k++)
        fmp[i][j] += ctf[j][k] * cmp[i][k];
    }
    for (j = 4; j < 10; j++) {
      fmp[i][j] = 0.0;
      /*for (k = 5; k < 10; k++)*/
      for (k = 4; k < 10; k++)
  fmp[i][j] += ctf[j][k] * cmp[i][k];
    }
  }
}

/* ----------------------------------------------------------------------
   cart_to_frac = Cartesian to fractional
   cart_to_frac computes a transformation matrix to convert
     a multipole object in Cartesian coordinates to fractional
   note the multipole components are stored in the condensed
     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
------------------------------------------------------------------------- */

void PairAmoeba::cart_to_frac()
{
  int i,j,k,m,i1,i2;
  int qi1[6] = {0,1,2,0,0,1};    // decremented vs Fortran
  int qi2[6] = {0,1,2,1,2,2};    // decremented vs Fortran
  double a[3][3];   // indices not flipped vs Fortran

  // set the reciprocal vector transformation matrix

  for (i = 0; i < 3; i++) {
    a[0][i] = nfft1 * recip[i][0];
    a[1][i] = nfft2 * recip[i][1];
    a[2][i] = nfft3 * recip[i][2];
  }

  // get the Cartesian to fractional conversion matrix

  for (i = 0; i < 10; i++)
    for (j = 0; j < 10; j++)
      ctf[i][j] = 0.0;

  ctf[0][0] = 1.0;
  for (i = 1; i < 4; i++)
    for (j = 1; j < 4; j++)
      ctf[i][j] = a[i-1][j-1];

  for (i1 = 0; i1 < 3; i1++) {
    k = qi1[i1];
    for (i2 = 0; i2 < 6; i2++) {
      i = qi1[i2];
      j = qi2[i2];
      ctf[i1+4][i2+4] = a[k][i] * a[k][j];
    }
  }

  for (i1 = 3; i1 < 6; i1++) {
    k = qi1[i1];
    m = qi2[i1];
    for (i2 = 0; i2 < 6; i2++) {
      i = qi1[i2];
      j = qi2[i2];
      ctf[i1+4][i2+4] = a[k][i]*a[m][j] + a[k][j]*a[m][i];
    }
  }
}

/* ----------------------------------------------------------------------
   fphi_to_cphi = transformation of potential
   fphi_to_cphi transforms the reciprocal space potential from
   fractional to Cartesian coordinates
------------------------------------------------------------------------- */

void PairAmoeba::fphi_to_cphi(double **fphi, double **cphi)
{
  int i,j,k;

  // find the matrix to convert fractional to Cartesian

  frac_to_cart();

  // apply the transformation to get the Cartesian potential

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    cphi[i][0] = ftc[0][0] * fphi[i][0];
    for (j = 1; j < 4; j++) {
      cphi[i][j] = 0.0;
      for (k = 1; k < 4; k++)
        cphi[i][j] += ftc[j][k] * fphi[i][k];
    }
    for (j = 4; j < 10; j++) {
      cphi[i][j] = 0.0;
      for (k = 4; k < 10; k++)
        cphi[i][j] += ftc[j][k] * fphi[i][k];
    }
  }
}

/* ----------------------------------------------------------------------
   frac_to_cart = fractional to Cartesian
   frac_to_cart computes a transformation matrix to convert
     a multipole object in fraction coordinates to Cartesian
   note the multipole components are stored in the condensed
     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
------------------------------------------------------------------------- */

void PairAmoeba::frac_to_cart()
{
  int i,j,k,m,i1,i2;
  int qi1[6] = {0,1,2,0,0,1};    // decremented vs Fortran
  int qi2[6] = {0,1,2,1,2,2};    // decremented vs Fortran
  double a[3][3];   // indices not flipped vs Fortran

  // set the reciprocal vector transformation matrix

  for (i = 0; i < 3; i++) {
    a[i][0] = nfft1 * recip[i][0];
    a[i][1] = nfft2 * recip[i][1];
    a[i][2] = nfft3 * recip[i][2];
  }

  // get the fractional to Cartesian conversion matrix

  for (i = 0; i < 10; i++)
    for (j = 0; j < 10; j++)
      ftc[i][j] = 0.0;

  ftc[0][0] = 1.0;
  for (i = 1; i < 4; i++)
    for (j = 1; j < 4; j++)
      ftc[i][j] = a[i-1][j-1];

  for (i1 = 0; i1 < 3; i1++) {
    k = qi1[i1];
    for (i2 = 0; i2 < 3; i2++) {
      i = qi1[i2];
      ftc[i1+4][i2+4] = a[k][i] * a[k][i];
    }
    for (i2 = 3; i2 < 6; i2++) {
      i = qi1[i2];
      j = qi2[i2];
      ftc[i1+4][i2+4] = 2.0 * a[k][i] * a[k][j];
    }
  }

  for (i1 = 3; i1 < 6; i1++) {
    k = qi1[i1];
    m = qi2[i1];
    for (i2 = 0; i2 < 3; i2++) {
      i = qi1[i2];
      ftc[i1+4][i2+4] = a[k][i] * a[m][i];
    }
    for (i2 = 3; i2 < 6; i2++) {
      i = qi1[i2];
      j = qi2[i2];
      ftc[i1+4][i2+4] = a[k][i]*a[m][j] + a[m][i]*a[k][j];
    }
  }
}

/* ----------------------------------------------------------------------
   grid_mpole = put multipoles on PME grid
   grid_mpole maps fractional atomic multipoles to PME grid
------------------------------------------------------------------------- */

void PairAmoeba::grid_mpole(double **fmp, FFT_SCALAR ***grid)
{
  int i,j,k,m,ib,jb,kb;
  double v0,u0,t0;
  double v1,u1,t1;
  double v2,u2,t2;
  double term0,term1,term2;

  int nlpts = (bsorder-1) / 2;

  // spread the permanent multipole moments onto the grid

  int nlocal = atom->nlocal;

  for (m = 0; m < nlocal; m++) {

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      v0 = thetai3[m][kb][0];
      v1 = thetai3[m][kb][1];
      v2 = thetai3[m][kb][2];

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        u0 = thetai2[m][jb][0];
        u1 = thetai2[m][jb][1];
        u2 = thetai2[m][jb][2];
        term0 = fmp[m][0]*u0*v0 + fmp[m][2]*u1*v0 + fmp[m][3]*u0*v1 +
          fmp[m][5]*u2*v0 + fmp[m][6]*u0*v2 + fmp[m][9]*u1*v1;
        term1 = fmp[m][1]*u0*v0 + fmp[m][7]*u1*v0 + fmp[m][8]*u0*v1;
        term2 = fmp[m][4]*u0*v0;

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          t0 = thetai1[m][ib][0];
          t1 = thetai1[m][ib][1];
          t2 = thetai1[m][ib][2];
          grid[k][j][i] += term0*t0 + term1*t1 + term2*t2;

          // if (m == 0) {
          //   int istencil = kb*bsorder*bsorder + jb*bsorder + ib + 1;
          //   printf("GRIDMPOLE iStencil %d atomID %d igrid %d %d %d "
                //    "ibjbkb %d %d %d ijk %d %d %d "
                //    "th1 %g %g %g th2 %g %g %g th3 %g %g %g "
                //    "fmp0-9 %e %e %e %e %e %e %e %e %e %e "
    //    "term012 %e %e %e "
                //    "gridvalue %g\n",
                //    istencil,atom->tag[m],
                //    igrid[m][0],igrid[m][1],igrid[m][2],
                //    ib,jb,kb,i,j,k,
                //    t0,t1,t2,u0,u1,u2,v0,v1,v2,
                //    fmp[m][0],fmp[m][1],fmp[m][2],
                //    fmp[m][3],fmp[m][4],fmp[m][5],
                //    fmp[m][6],fmp[m][7],fmp[m][8],fmp[m][9],
    //    term0,term1,term2,
                //    term0*t0 + term1*t1 + term2*t2);
          // }

          i++;
        }
        j++;
      }
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   fphi_mpole = multipole potential from grid
   fphi_mpole extracts the permanent multipole potential from
   the particle mesh Ewald grid
------------------------------------------------------------------------- */

void PairAmoeba::fphi_mpole(FFT_SCALAR ***grid, double **fphi)
{
  int i,j,k,m,ib,jb,kb;
  double v0,v1,v2,v3;
  double u0,u1,u2,u3;
  double t0,t1,t2,t3,tq;
  double tu00,tu10,tu01,tu20,tu11;
  double tu02,tu21,tu12,tu30,tu03;
  double tuv000,tuv100,tuv010,tuv001;
  double tuv200,tuv020,tuv002,tuv110;
  double tuv101,tuv011,tuv300,tuv030;
  double tuv003,tuv210,tuv201,tuv120;
  double tuv021,tuv102,tuv012,tuv111;

  int nlpts = (bsorder-1) / 2;

  // extract the permanent multipole field at each site

  int nlocal = atom->nlocal;

  for (m = 0; m < nlocal; m++) {
    tuv000 = 0.0;
    tuv001 = 0.0;
    tuv010 = 0.0;
    tuv100 = 0.0;
    tuv200 = 0.0;
    tuv020 = 0.0;
    tuv002 = 0.0;
    tuv110 = 0.0;
    tuv101 = 0.0;
    tuv011 = 0.0;
    tuv300 = 0.0;
    tuv030 = 0.0;
    tuv003 = 0.0;
    tuv210 = 0.0;
    tuv201 = 0.0;
    tuv120 = 0.0;
    tuv021 = 0.0;
    tuv102 = 0.0;
    tuv012 = 0.0;
    tuv111 = 0.0;

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      v0 = thetai3[m][kb][0];
      v1 = thetai3[m][kb][1];
      v2 = thetai3[m][kb][2];
      v3 = thetai3[m][kb][3];
      tu00 = 0.0;
      tu10 = 0.0;
      tu01 = 0.0;
      tu20 = 0.0;
      tu11 = 0.0;
      tu02 = 0.0;
      tu30 = 0.0;
      tu21 = 0.0;
      tu12 = 0.0;
      tu03 = 0.0;

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        u0 = thetai2[m][jb][0];
        u1 = thetai2[m][jb][1];
        u2 = thetai2[m][jb][2];
        u3 = thetai2[m][jb][3];
        t0 = 0.0;
        t1 = 0.0;
        t2 = 0.0;
        t3 = 0.0;

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          tq = grid[k][j][i];
          t0 += tq*thetai1[m][ib][0];
          t1 += tq*thetai1[m][ib][1];
          t2 += tq*thetai1[m][ib][2];
          t3 += tq*thetai1[m][ib][3];
          i++;
        }

        tu00 += t0*u0;
        tu10 += t1*u0;
        tu01 += t0*u1;
        tu20 += t2*u0;
        tu11 += t1*u1;
        tu02 += t0*u2;
        tu30 += t3*u0;
        tu21 += t2*u1;
        tu12 += t1*u2;
        tu03 += t0*u3;
        j++;
      }

      tuv000 += tu00*v0;
      tuv100 += tu10*v0;
      tuv010 += tu01*v0;
      tuv001 += tu00*v1;
      tuv200 += tu20*v0;
      tuv020 += tu02*v0;
      tuv002 += tu00*v2;
      tuv110 += tu11*v0;
      tuv101 += tu10*v1;
      tuv011 += tu01*v1;
      tuv300 += tu30*v0;
      tuv030 += tu03*v0;
      tuv003 += tu00*v3;
      tuv210 += tu21*v0;
      tuv201 += tu20*v1;
      tuv120 += tu12*v0;
      tuv021 += tu02*v1;
      tuv102 += tu10*v2;
      tuv012 += tu01*v2;
      tuv111 += tu11*v1;
      k++;
    }

    fphi[m][0] = tuv000;
    fphi[m][1] = tuv100;
    fphi[m][2] = tuv010;
    fphi[m][3] = tuv001;
    fphi[m][4] = tuv200;
    fphi[m][5] = tuv020;
    fphi[m][6] = tuv002;
    fphi[m][7] = tuv110;
    fphi[m][8] = tuv101;
    fphi[m][9] = tuv011;
    fphi[m][10] = tuv300;
    fphi[m][11] = tuv030;
    fphi[m][12] = tuv003;
    fphi[m][13] = tuv210;
    fphi[m][14] = tuv201;
    fphi[m][15] = tuv120;
    fphi[m][16] = tuv021;
    fphi[m][17] = tuv102;
    fphi[m][18] = tuv012;
    fphi[m][19] = tuv111;
  }
}

/* ----------------------------------------------------------------------
   grid_uind = put induced dipoles on PME grid
   grid_uind maps fractional induced dipoles to the PME grid
------------------------------------------------------------------------- */

void PairAmoeba::grid_uind(double **fuind, double **fuinp, FFT_SCALAR ****grid)
{
  int i,j,k,m,ib,jb,kb;
  double v0,u0,t0;
  double v1,u1,t1;
  double term01,term11;
  double term02,term12;

  int nlpts = (bsorder-1) / 2;

  // put the induced dipole moments onto the grid

  int nlocal = atom->nlocal;

  for (m = 0; m < nlocal; m++) {

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      v0 = thetai3[m][kb][0];
      v1 = thetai3[m][kb][1];

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        u0 = thetai2[m][jb][0];
        u1 = thetai2[m][jb][1];
        term01 = fuind[m][1]*u1*v0 + fuind[m][2]*u0*v1;
        term11 = fuind[m][0]*u0*v0;
        term02 = fuinp[m][1]*u1*v0 + fuinp[m][2]*u0*v1;
        term12 = fuinp[m][0]*u0*v0;

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          t0 = thetai1[m][ib][0];
          t1 = thetai1[m][ib][1];
          grid[k][j][i][0] += term01*t0 + term11*t1;
          grid[k][j][i][1] += term02*t0 + term12*t1;
    //printf("gridcheck %g %g \n",grid[k][j][i][0],grid[k][j][i][0]);
          i++;
        }
        j++;
      }
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   fphi_uind = induced potential from grid
   fphi_uind extracts the induced dipole potential from the particle mesh Ewald grid
------------------------------------------------------------------------- */

void PairAmoeba::fphi_uind(FFT_SCALAR ****grid, double **fdip_phi1,
                           double **fdip_phi2, double **fdip_sum_phi)
{
  int i,j,k,m,ib,jb,kb;
  double v0,v1,v2,v3;
  double u0,u1,u2,u3;
  double t0,t1,t2,t3;
  double t0_1,t0_2,t1_1,t1_2;
  double t2_1,t2_2,tq_1,tq_2;
  double tu00,tu10,tu01,tu20,tu11;
  double tu02,tu30,tu21,tu12,tu03;
  double tu00_1,tu01_1,tu10_1;
  double tu00_2,tu01_2,tu10_2;
  double tu20_1,tu11_1,tu02_1;
  double tu20_2,tu11_2,tu02_2;
  double tuv100_1,tuv010_1,tuv001_1;
  double tuv100_2,tuv010_2,tuv001_2;
  double tuv200_1,tuv020_1,tuv002_1;
  double tuv110_1,tuv101_1,tuv011_1;
  double tuv200_2,tuv020_2,tuv002_2;
  double tuv110_2,tuv101_2,tuv011_2;
  double tuv000,tuv100,tuv010,tuv001;
  double tuv200,tuv020,tuv002,tuv110;
  double tuv101,tuv011,tuv300,tuv030;
  double tuv003,tuv210,tuv201,tuv120;
  double tuv021,tuv102,tuv012,tuv111;

  int nlpts = (bsorder-1) / 2;

  // extract the permanent multipole field at each site

  int nlocal = atom->nlocal;

  for (m = 0; m < nlocal; m++) {
    tuv100_1 = 0.0;
    tuv010_1 = 0.0;
    tuv001_1 = 0.0;
    tuv200_1 = 0.0;
    tuv020_1 = 0.0;
    tuv002_1 = 0.0;
    tuv110_1 = 0.0;
    tuv101_1 = 0.0;
    tuv011_1 = 0.0;
    tuv100_2 = 0.0;
    tuv010_2 = 0.0;
    tuv001_2 = 0.0;
    tuv200_2 = 0.0;
    tuv020_2 = 0.0;
    tuv002_2 = 0.0;
    tuv110_2 = 0.0;
    tuv101_2 = 0.0;
    tuv011_2 = 0.0;
    tuv000 = 0.0;
    tuv001 = 0.0;
    tuv010 = 0.0;
    tuv100 = 0.0;
    tuv200 = 0.0;
    tuv020 = 0.0;
    tuv002 = 0.0;
    tuv110 = 0.0;
    tuv101 = 0.0;
    tuv011 = 0.0;
    tuv300 = 0.0;
    tuv030 = 0.0;
    tuv003 = 0.0;
    tuv210 = 0.0;
    tuv201 = 0.0;
    tuv120 = 0.0;
    tuv021 = 0.0;
    tuv102 = 0.0;
    tuv012 = 0.0;
    tuv111 = 0.0;

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      v0 = thetai3[m][kb][0];
      v1 = thetai3[m][kb][1];
      v2 = thetai3[m][kb][2];
      v3 = thetai3[m][kb][3];
      tu00_1 = 0.0;
      tu01_1 = 0.0;
      tu10_1 = 0.0;
      tu20_1 = 0.0;
      tu11_1 = 0.0;
      tu02_1 = 0.0;
      tu00_2 = 0.0;
      tu01_2 = 0.0;
      tu10_2 = 0.0;
      tu20_2 = 0.0;
      tu11_2 = 0.0;
      tu02_2 = 0.0;
      tu00 = 0.0;
      tu10 = 0.0;
      tu01 = 0.0;
      tu20 = 0.0;
      tu11 = 0.0;
      tu02 = 0.0;
      tu30 = 0.0;
      tu21 = 0.0;
      tu12 = 0.0;
      tu03 = 0.0;

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        u0 = thetai2[m][jb][0];
        u1 = thetai2[m][jb][1];
        u2 = thetai2[m][jb][2];
        u3 = thetai2[m][jb][3];
        t0_1 = 0.0;
        t1_1 = 0.0;
        t2_1 = 0.0;
        t0_2 = 0.0;
        t1_2 = 0.0;
        t2_2 = 0.0;
        t3 = 0.0;

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          tq_1 = grid[k][j][i][0];
          tq_2 = grid[k][j][i][1];
          t0_1 += tq_1*thetai1[m][ib][0];
          t1_1 += tq_1*thetai1[m][ib][1];
          t2_1 += tq_1*thetai1[m][ib][2];
          t0_2 += tq_2*thetai1[m][ib][0];
          t1_2 += tq_2*thetai1[m][ib][1];
          t2_2 += tq_2*thetai1[m][ib][2];
          t3 += (tq_1+tq_2)*thetai1[m][ib][3];
          i++;
        }

        tu00_1 += t0_1*u0;
        tu10_1 += t1_1*u0;
        tu01_1 += t0_1*u1;
        tu20_1 += t2_1*u0;
        tu11_1 += t1_1*u1;
        tu02_1 += t0_1*u2;
        tu00_2 += t0_2*u0;
        tu10_2 += t1_2*u0;
        tu01_2 += t0_2*u1;
        tu20_2 += t2_2*u0;
        tu11_2 += t1_2*u1;
        tu02_2 += t0_2*u2;
        t0 = t0_1 + t0_2;
        t1 = t1_1 + t1_2;
        t2 = t2_1 + t2_2;
        tu00 += t0*u0;
        tu10 += t1*u0;
        tu01 += t0*u1;
        tu20 += t2*u0;
        tu11 += t1*u1;
        tu02 += t0*u2;
        tu30 += t3*u0;
        tu21 += t2*u1;
        tu12 += t1*u2;
        tu03 += t0*u3;
        j++;
      }

      tuv100_1 += tu10_1*v0;
      tuv010_1 += tu01_1*v0;
      tuv001_1 += tu00_1*v1;
      tuv200_1 += tu20_1*v0;
      tuv020_1 += tu02_1*v0;
      tuv002_1 += tu00_1*v2;
      tuv110_1 += tu11_1*v0;
      tuv101_1 += tu10_1*v1;
      tuv011_1 += tu01_1*v1;
      tuv100_2 += tu10_2*v0;
      tuv010_2 += tu01_2*v0;
      tuv001_2 += tu00_2*v1;
      tuv200_2 += tu20_2*v0;
      tuv020_2 += tu02_2*v0;
      tuv002_2 += tu00_2*v2;
      tuv110_2 += tu11_2*v0;
      tuv101_2 += tu10_2*v1;
      tuv011_2 += tu01_2*v1;
      tuv000 += tu00*v0;
      tuv100 += tu10*v0;
      tuv010 += tu01*v0;
      tuv001 += tu00*v1;
      tuv200 += tu20*v0;
      tuv020 += tu02*v0;
      tuv002 += tu00*v2;
      tuv110 += tu11*v0;
      tuv101 += tu10*v1;
      tuv011 += tu01*v1;
      tuv300 += tu30*v0;
      tuv030 += tu03*v0;
      tuv003 += tu00*v3;
      tuv210 += tu21*v0;
      tuv201 += tu20*v1;
      tuv120 += tu12*v0;
      tuv021 += tu02*v1;
      tuv102 += tu10*v2;
      tuv012 += tu01*v2;
      tuv111 += tu11*v1;
      k++;
    }

    fdip_phi1[m][0] = 0.0;
    fdip_phi1[m][1] = tuv100_1;
    fdip_phi1[m][2] = tuv010_1;
    fdip_phi1[m][3] = tuv001_1;
    fdip_phi1[m][4] = tuv200_1;
    fdip_phi1[m][5] = tuv020_1;
    fdip_phi1[m][6] = tuv002_1;
    fdip_phi1[m][7] = tuv110_1;
    fdip_phi1[m][8] = tuv101_1;
    fdip_phi1[m][9] = tuv011_1;

    fdip_phi2[m][0] = 0.0;
    fdip_phi2[m][1] = tuv100_2;
    fdip_phi2[m][2] = tuv010_2;
    fdip_phi2[m][3] = tuv001_2;
    fdip_phi2[m][4] = tuv200_2;
    fdip_phi2[m][5] = tuv020_2;
    fdip_phi2[m][6] = tuv002_2;
    fdip_phi2[m][7] = tuv110_2;
    fdip_phi2[m][8] = tuv101_2;
    fdip_phi2[m][9] = tuv011_2;

    fdip_sum_phi[m][0] = tuv000;
    fdip_sum_phi[m][1] = tuv100;
    fdip_sum_phi[m][2] = tuv010;
    fdip_sum_phi[m][3] = tuv001;
    fdip_sum_phi[m][4] = tuv200;
    fdip_sum_phi[m][5] = tuv020;
    fdip_sum_phi[m][6] = tuv002;
    fdip_sum_phi[m][7] = tuv110;
    fdip_sum_phi[m][8] = tuv101;
    fdip_sum_phi[m][9] = tuv011;
    fdip_sum_phi[m][10] = tuv300;
    fdip_sum_phi[m][11] = tuv030;
    fdip_sum_phi[m][12] = tuv003;
    fdip_sum_phi[m][13] = tuv210;
    fdip_sum_phi[m][14] = tuv201;
    fdip_sum_phi[m][15] = tuv120;
    fdip_sum_phi[m][16] = tuv021;
    fdip_sum_phi[m][17] = tuv102;
    fdip_sum_phi[m][18] = tuv012;
    fdip_sum_phi[m][19] = tuv111;
  }
}

/* ----------------------------------------------------------------------
   grid_disp = put dispersion sites on PME grid
   grid_disp maps dispersion coefficients to PME grid
------------------------------------------------------------------------- */

void PairAmoeba::grid_disp(FFT_SCALAR ***grid)
{
  int i,j,k,m,ib,jb,kb,itype,iclass;
  double v0,u0,t0;
  double term;

  int nlpts = (bsorder-1) / 2;

  // put the dispersion sites onto the grid

  int nlocal = atom->nlocal;

  for (m = 0; m < nlocal; m++) {
    itype = amtype[m];
    iclass = amtype2class[itype];

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      v0 = thetai3[m][kb][0] * csix[iclass];

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        u0 = thetai2[m][jb][0];
        term = v0 * u0;

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          t0 = thetai1[m][ib][0];
          grid[k][j][i] += term*t0;
          i++;
        }
        j++;
      }
      k++;
    }
  }
}

/* ----------------------------------------------------------------------
   kewald = setup for particle mesh Ewald sum
   kewald assigns particle mesh Ewald parameters and options for a periodic system
------------------------------------------------------------------------- */

void PairAmoeba::kewald()
{
  int nfft1,nfft2,nfft3;
  double delta;
  double edens,ddens;
  double size,slope;

  // use_ewald is for multipole and induce and polar

  if (!use_ewald) aeewald = apewald = 0.0;

  if (use_ewald) {
    if (!aeewald_key) aeewald = ewaldcof(ewaldcut);

    if (!apewald_key) {
      apewald = aeewald;
      size = MIN(domain->xprd,domain->yprd);
      size = MIN(size,domain->zprd);
      if (size < 6.0) {
        slope = (1.0-apewald) / 2.0;
        apewald += slope*(6.0-size);
      }
    }

    if (!pmegrid_key) {
      delta = 1.0e-8;
      edens = 1.2;
      nefft1 = static_cast<int> (domain->xprd*edens-delta) + 1;
      nefft2 = static_cast<int> (domain->yprd*edens-delta) + 1;
      nefft3 = static_cast<int> (domain->zprd*edens-delta) + 1;
    }

    // round grid sizes up to be factorable by 2,3,5
    // NOTE: also worry about satisfying Tinker minfft ?

    while (!factorable(nefft1)) nefft1++;
    while (!factorable(nefft2)) nefft2++;
    while (!factorable(nefft3)) nefft3++;
  }

  // use_dewald is for dispersion

  if (!use_dewald) adewald = 0.0;

  if (use_dewald) {
    if (!adewald_key) adewald = ewaldcof(dispcut);

    if (!dpmegrid_key) {
      delta = 1.0e-8;
      ddens = 0.8;
      ndfft1 = static_cast<int> (domain->xprd*ddens-delta) + 1;
      ndfft2 = static_cast<int> (domain->yprd*ddens-delta) + 1;
      ndfft3 = static_cast<int> (domain->zprd*ddens-delta) + 1;
    }

    // increase grid sizes until factorable by 2,3,5
    // NOTE: also worry about satisfying Tinker minfft ?

    while (!factorable(ndfft1)) ndfft1++;
    while (!factorable(ndfft2)) ndfft2++;
    while (!factorable(ndfft3)) ndfft3++;
  }

  // done if no Ewald

  if (!use_ewald && !use_dewald) return;

  // set maximum sizes for PME grids and B-spline order
  // allocate vectors and arrays accordingly

  nfft1 = nfft2 = nfft3 = 0;
  if (use_ewald) nfft1 = nefft1;
  if (use_ewald) nfft2 = nefft2;
  if (use_ewald) nfft3 = nefft3;
  if (use_dewald) nfft1 = MAX(nfft1,ndfft1);
  if (use_dewald) nfft2 = MAX(nfft2,ndfft2);
  if (use_dewald) nfft3 = MAX(nfft3,ndfft3);

  bsordermax = 0;
  if (use_ewald) bsordermax = bseorder;
  if (use_ewald) bsordermax = MAX(bsordermax,bsporder);
  if (use_dewald) bsordermax = MAX(bsordermax,bsdorder);

  memory->create(bsmod1,nfft1,"amoeba:bsmod1");
  memory->create(bsmod2,nfft2,"amoeba:bsmod2");
  memory->create(bsmod3,nfft3,"amoeba:bsmod3");
  memory->create(bsbuild,bsordermax,bsordermax,"amoeba:bsbuild");
}

/* ----------------------------------------------------------------------
   ewaldcof = estimation of Ewald coefficient
   ewaldcof finds an Ewald coefficient such that all terms
   beyond the specified cutoff distance will have a value less
   than a specified tolerance
------------------------------------------------------------------------- */

double PairAmoeba::ewaldcof(double cutoff)
{
  int i,k,m;
  double x,xlo,xhi,y;
  double eps,ratio;

  // set the tolerance value; use of 1.0d-8 instead of 1.0d-6
  // gives large coefficients that ensure gradient continuity

  eps = 1.0e-8;

  // get approximate value from cutoff and tolerance

  ratio = eps + 1.0;
  x = 0.5;
  i = 0;
  while (ratio >= eps) {
    i++;
    x *= 2.0;
    y = x * cutoff;
    ratio = erfc(y) / cutoff;
  }

  // use a binary search to refine the coefficient

  k = i + 60;
  xlo = 0.0;
  xhi = x;
  for (m = 0; m < k; m++) {
    x = (xlo+xhi) / 2.0;
    y = x * cutoff;
    ratio = erfc(y) / cutoff;
    if (ratio >= eps) xlo = x;
    else xhi = x;
  }

  return x;
}

/* ----------------------------------------------------------------------
   check if all factors of n are in factors[] list of length nfactors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int PairAmoeba::factorable(int n)
{
  int i;

  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        break;
      }
    }
    if (i == nfactors) return 0;
  }

  return 1;
}
