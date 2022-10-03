// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_amoeba.h"

#include "amoeba_convolution.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

using MathSpecial::cube;
using MathSpecial::powint;

enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   dispersion = Ewald dispersion
   adapted from Tinker edisp1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::dispersion()
{
  // set cutoffs, taper coeffs, and PME params

  if (use_dewald) choose(DISP_LONG);
  else choose(DISP);

  // owned atoms

  int nlocal = atom->nlocal;

  // compute the real space portion of the Ewald summation

  if (disp_rspace_flag) dispersion_real();

  // compute the reciprocal space part of the Ewald summation

  if (disp_kspace_flag) dispersion_kspace();

  // compute the self-energy portion of the Ewald summation

  int itype,iclass;
  double term;

  for (int i = 0; i < nlocal; i++) {
    itype = amtype[i];
    iclass = amtype2class[itype];
    term = powint(aewald,6) / 12.0;
    edisp += term*csix[iclass]*csix[iclass];
  }
}

/* ----------------------------------------------------------------------
   dispersion_real = real-space portion of Ewald dispersion
   adapted from Tinker edreal1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::dispersion_real()
{
  int i,j,ii,jj,itype,jtype,iclass,jclass;
  double xi,yi,zi;
  double xr,yr,zr;
  double e,de;
  double ci,ck;
  double r,r2,r6,r7;
  double ai,ai2;
  double ak,ak2;
  double di,di2,di3,di4,di5;
  double dk,dk2,dk3;
  double ti,ti2;
  double tk,tk2;
  double expi,expk;
  double damp3,damp5;
  double damp,ddamp;
  double ralpha2,scale;
  double expterm,term;
  double expa,rterm;
  double dedx,dedy,dedz;
  double vxx,vyx,vzx;
  double vyy,vzy,vzz;
  double factor_disp;

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute the real space portion of the Ewald summation

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    ci = csix[iclass];
    ai = adisp[iclass];
    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];

    // decide whether to compute the current interaction

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_disp = special_disp[sbmask15(j)];
      j &= NEIGHMASK15;

      xr = xi - x[j][0];
      yr = yi - x[j][1];
      zr = zi - x[j][2];
      r2 = xr*xr + yr*yr + zr*zr;
      if (r2 > off2) continue;

      // compute the energy contribution for this interaction

      jtype = amtype[j];
      jclass = amtype2class[jtype];
      ck = csix[jclass];
      ak = adisp[jclass];

      r6 = r2*r2*r2;
      ralpha2 = r2 * aewald*aewald;
      term = 1.0 + ralpha2 + 0.5*ralpha2*ralpha2;
      expterm = exp(-ralpha2);
      expa = expterm * term;

      // find the damping factor for the dispersion interaction

      r = sqrt(r2);
      r7 = r6 * r;
      di = ai * r;
      di2 = di * di;
      di3 = di * di2;
      dk = ak * r;
      expi = exp(-di);
      expk = exp(-dk);

      if (ai != ak) {
        ai2 = ai * ai;
        ak2 = ak * ak;
        dk2 = dk * dk;
        dk3 = dk * dk2;
        ti = ak2 / (ak2-ai2);
        ti2 = ti * ti;
        tk = ai2 / (ai2-ak2);
        tk2 = tk * tk;
        damp3 = 1.0 - ti2*(1.0+di+0.5*di2)*expi - tk2*(1.0+dk+0.5*dk2)*expk -
          2.0*ti2*tk*(1.0+di)*expi - 2.0*tk2*ti*(1.0+dk)*expk;
        damp5 = 1.0 - ti2*(1.0+di+0.5*di2+di3/6.0)*expi -
          tk2*(1.0+dk+0.5*dk2 + dk3/6.0)*expk -
          2.0*ti2*tk*(1.0+di+di2/3.0)*expi - 2.0*tk2*ti*(1.0+dk+dk2/3.0)*expk;
        ddamp = 0.25 * di2 * ti2 * ai * expi * (r*ai+4.0*tk-1.0) +
          0.25 * dk2 * tk2 * ak * expk * (r*ak+4.0*ti-1.0);

      } else {
        di4 = di2 * di2;
        di5 = di2 * di3;
        damp3 = 1.0 - (1.0+di+0.5*di2 + 7.0*di3/48.0+di4/48.0)*expi;
        damp5 = 1.0 - (1.0+di+0.5*di2 + di3/6.0+di4/24.0+di5/144.0)*expi;
        ddamp = ai * expi * (di5-3.0*di3-3.0*di2) / 96.0;
      }

      damp = 1.5*damp5 - 0.5*damp3;

      // apply damping and scaling factors for this interaction

      scale = factor_disp * damp*damp;
      scale = scale - 1.0;
      e = -ci * ck * (expa+scale) / r6;
      rterm = -cube(ralpha2) * expterm / r;
      de = -6.0*e/r2 - ci*ck*rterm/r7 - 2.0*ci*ck*factor_disp*damp*ddamp/r7;

      edisp += e;

      // increment the damped dispersion derivative components

      dedx = de * xr;
      dedy = de * yr;
      dedz = de * zr;
      f[i][0] -= dedx;
      f[i][1] -= dedy;
      f[i][2] -= dedz;
      f[j][0] += dedx;
      f[j][1] += dedy;
      f[j][2] += dedz;

      // increment the internal virial tensor components

      if (vflag_global) {
        vxx = xr * dedx;
        vyx = yr * dedx;
        vzx = zr * dedx;
        vyy = yr * dedy;
        vzy = zr * dedy;
        vzz = zr * dedz;

        virdisp[0] -= vxx;
        virdisp[1] -= vyy;
        virdisp[2] -= vzz;
        virdisp[3] -= vyx;
        virdisp[4] -= vzx;
        virdisp[5] -= vzy;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   dispersion_kspace = KSpace portion of Ewald dispersion
   adapted from Tinker edrecip1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::dispersion_kspace()
{
  int i,j,k,m,n,ib,jb,kb,itype,iclass;
  int nhalf1,nhalf2,nhalf3;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  double e,fi,denom,scale;
  double r1,r2,r3;
  double h1,h2,h3;
  double term,vterm;
  double expterm;
  double erfcterm;
  double hsq,struc2;
  double h,hhh,b,bfac;
  double term1,denom0;
  double fac1,fac2,fac3;
  double de1,de2,de3;
  double dt1,dt2,dt3;
  double t1,t2,t3;

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  // owned atoms

  double **f = atom->f;
  int nlocal = atom->nlocal;

  double volbox = domain->prd[0] * domain->prd[1] * domain->prd[2];

  // FFT moduli pre-computations
  // set igrid for each atom and its B-spline coeffs

  nfft1 = d_kspace->nx;
  nfft2 = d_kspace->ny;
  nfft3 = d_kspace->nz;
  bsorder = d_kspace->order;

  moduli();
  bspline_fill();

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by zero()

  double ***gridpre = (double ***) d_kspace->zero();

  // map atoms to grid

  grid_disp(gridpre);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomposition

  double *gridfft = d_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  nhalf1 = (nfft1+1) / 2;
  nhalf2 = (nfft2+1) / 2;
  nhalf3 = (nfft3+1) / 2;

  nxlo = d_kspace->nxlo_fft;
  nxhi = d_kspace->nxhi_fft;
  nylo = d_kspace->nylo_fft;
  nyhi = d_kspace->nyhi_fft;
  nzlo = d_kspace->nzlo_fft;
  nzhi = d_kspace->nzhi_fft;

  bfac = MY_PI / aewald;
  fac1 = 2.0*pow(MY_PI,3.5);
  fac2 = cube(aewald);
  fac3 = -2.0*aewald*MY_PI*MY_PI;
  denom0 = (6.0*volbox)/pow(MY_PI,1.5);

  n = 0;
  for (k = nzlo; k <= nzhi; k++) {
    for (j = nylo; j <= nyhi; j++) {
      for (i = nxlo; i <= nxhi; i++) {
        r1 = (i >= nhalf1) ? i-nfft1 : i;
        r2 = (j >= nhalf2) ? j-nfft2 : j;
        r3 = (k >= nhalf3) ? k-nfft3 : k;
        h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;  // matvec
        h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
        h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
        hsq = h1*h1 + h2*h2 + h3*h3;
        h = sqrt(hsq);
        b = h*bfac;
        hhh = h*hsq;
        term = -b*b;
        expterm = 0.0;
        erfcterm = erfc(b);
        denom = denom0*bsmod1[i]*bsmod2[j]*bsmod3[k];
        if (term > -50.0 && hsq != 0.0) {
          expterm = exp(term);
          erfcterm = erfc(b);
          term1 = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq);
          struc2 = gridfft[n]*gridfft[n] + gridfft[n+1]*gridfft[n+1];
          e = -(term1 / denom) * struc2;
          edisp += e;
          if (vflag_global) {
            vterm = 3.0 * (fac1*erfcterm*h + fac3*expterm) * struc2/denom;
            virdisp[0] -= h1*h1*vterm - e;
            virdisp[1] -= h2*h2*vterm - e;
            virdisp[2] -= h3*h3*vterm - e;
            virdisp[3] -= h1*h2*vterm;
            virdisp[4] -= h1*h3*vterm;
            virdisp[5] -= h2*h3*vterm;
          }
        } else term1 = 0.0;
        scale = -term1 / denom;
        gridfft[n] *= scale;
        gridfft[n+1] *= scale;
        n += 2;
      }
    }
  }

  // post-convolution operations including backward FFT
  // gridppost = my portion of 3d grid in brick decomp w/ ghost values

  double ***gridpost = (double ***) d_kspace->post_convolution();

  // get first derivatives of the reciprocal space energy

  int nlpts = (bsorder-1) / 2;

  for (m = 0; m < nlocal; m++) {
    itype = amtype[m];
    iclass = amtype2class[itype];
    de1 = de2 = de3 = 0.0;

    k = igrid[m][2] - nlpts;
    for (kb = 0; kb < bsorder; kb++) {
      t3 = thetai3[m][kb][0];
      dt3 = nfft3 * thetai3[m][kb][1];

      j = igrid[m][1] - nlpts;
      for (jb = 0; jb < bsorder; jb++) {
        t2 = thetai2[m][jb][0];
        dt2 = nfft2 * thetai2[m][jb][1];

        i = igrid[m][0] - nlpts;
        for (ib = 0; ib < bsorder; ib++) {
          t1 = thetai1[m][ib][0];
          dt1 = nfft1 * thetai1[m][ib][1];
          term = gridpost[k][j][i];
          de1 += 2.0*term*dt1*t2*t3;
          de2 += 2.0*term*dt2*t1*t3;
          de3 += 2.0*term*dt3*t1*t2;
          i++;
        }
        j++;
      }
      k++;
    }

    fi = csix[iclass];
    f[m][0] -= fi * (recip[0][0]*de1 + recip[0][1]*de2 + recip[0][2]*de3);
    f[m][1] -= fi * (recip[1][0]*de1 + recip[1][1]*de2 + recip[1][2]*de3);
    f[m][2] -= fi * (recip[2][0]*de1 + recip[2][1]*de2 + recip[2][2]*de3);
  }

  // account for the energy and virial correction terms

  term = csixpr * aewald*aewald*aewald / denom0;

  if (comm->me == 0) {
    edisp -= term;
    if (vflag_global) {
      virdisp[0] -= term;
      virdisp[1] -= term;
      virdisp[2] -= term;
    }
  }
}
