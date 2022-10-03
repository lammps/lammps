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
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

using MathSpecial::square;

enum{FIELD,ZRSD,TORQUE,UFLD};                          // reverse comm
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ----------------------------------------------------------------------
   multipole = multipole interactions
   adapted from Tinker empole1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::multipole()
{
  double e;
  double felec;
  double term,fterm;
  double ci;
  double dix,diy,diz;
  double qixx,qixy,qixz,qiyy,qiyz,qizz;
  double cii,dii,qii;

  // set cutoffs, taper coeffs, and PME params

  if (use_ewald) choose(MPOLE_LONG);
  else choose(MPOLE);

  // owned atoms

  const int nlocal = atom->nlocal;

  // zero repulsion torque on owned + ghost atoms

  const int nall = nlocal + atom->nghost;

  for (int i = 0; i < nall; i++) {
    tq[i][0] = 0.0;
    tq[i][1] = 0.0;
    tq[i][2] = 0.0;
  }

  // set the energy unit conversion factor

  felec = electric / am_dielectric;

  // compute the real space part of the Ewald summation

  if (mpole_rspace_flag) multipole_real();

  // compute the reciprocal space part of the Ewald summation

  if (mpole_kspace_flag) multipole_kspace();

  // compute the Ewald self-energy term over all the atoms

  term = 2.0 * aewald * aewald;
  fterm = -felec * aewald / MY_PIS;

  for (int i = 0; i < nlocal; i++) {
    ci = rpole[i][0];
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];
    cii = ci*ci;
    dii = dix*dix + diy*diy + diz*diz;
    qii = 2.0*(qixy*qixy+qixz*qixz+qiyz*qiyz) +
      qixx*qixx + qiyy*qiyy + qizz*qizz;
    e = fterm * (cii + term*(dii/3.0+2.0*term*qii/5.0));
    empole += e;
  }
}

/* ----------------------------------------------------------------------
   multipole_real = real-space portion of mulipole interactions
   adapted from Tinker emreal1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::multipole_real()
{
  int i,j,k,itype,jtype,iclass,jclass;
  int ii,jj;
  int ix,iy,iz;
  double e,de,felec;
  double bfac;
  double alsq2,alsq2n;
  double exp2a,ralpha;
  double scalek;
  double xi,yi,zi;
  double xr,yr,zr;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double r,r2,rr1,rr3;
  double rr5,rr7,rr9,rr11;
  double rr1i,rr3i,rr5i,rr7i;
  double rr1k,rr3k,rr5k,rr7k;
  double rr1ik,rr3ik,rr5ik;
  double rr7ik,rr9ik,rr11ik;
  double ci,dix,diy,diz;
  double qixx,qixy,qixz;
  double qiyy,qiyz,qizz;
  double ck,dkx,dky,dkz;
  double qkxx,qkxy,qkxz;
  double qkyy,qkyz,qkzz;
  double dir,dkr,dik,qik;
  double qix,qiy,qiz,qir;
  double qkx,qky,qkz,qkr;
  double diqk,dkqi,qiqk;
  double dirx,diry,dirz;
  double dkrx,dkry,dkrz;
  double dikx,diky,dikz;
  double qirx,qiry,qirz;
  double qkrx,qkry,qkrz;
  double qikx,qiky,qikz;
  double qixk,qiyk,qizk;
  double qkxi,qkyi,qkzi;
  double qikrx,qikry,qikrz;
  double qkirx,qkiry,qkirz;
  double diqkx,diqky,diqkz;
  double dkqix,dkqiy,dkqiz;
  double diqkrx,diqkry,diqkrz;
  double dkqirx,dkqiry,dkqirz;
  double dqikx,dqiky,dqikz;
  double corei,corek;
  double vali,valk;
  double alphai,alphak;
  double term1,term2,term3;
  double term4,term5,term6;
  double term1i,term2i,term3i;
  double term1k,term2k,term3k;
  double term1ik,term2ik,term3ik;
  double term4ik,term5ik;
  double frcx,frcy,frcz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double factor_mpole;
  double ttmi[3],ttmk[3];
  double fix[3],fiy[3],fiz[3];
  double dmpi[9],dmpj[9];
  double dmpij[11];
  double bn[6];

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // owned atoms

  double *pval = atom->dvector[index_pval];
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // set conversion factor, cutoff and switching coefficients

  felec = electric / am_dielectric;

  // DEBUG

  //int count = 0;
  //int imin,imax;

  // compute the real space portion of the Ewald summation

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];
    ci = rpole[i][0];
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];
    if (!amoeba) {
      corei = pcore[iclass];
      alphai = palpha[iclass];
      vali = pval[i];
    }

    // evaluate all sites within the cutoff distance

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mpole = special_mpole[sbmask15(j)];
      j &= NEIGHMASK15;

      xr = x[j][0] - xi;
      yr = x[j][1] - yi;
      zr = x[j][2] - zi;
      r2 = xr*xr + yr*yr + zr*zr;
      if (r2 > off2) continue;

      // DEBUG

      //imin = MIN(atom->tag[i],atom->tag[j]);
      //imax = MAX(atom->tag[i],atom->tag[j]);

      jtype = amtype[j];
      jclass = amtype2class[jtype];

      r = sqrt(r2);
      ck = rpole[j][0];
      dkx = rpole[j][1];
      dky = rpole[j][2];
      dkz = rpole[j][3];
      qkxx = rpole[j][4];
      qkxy = rpole[j][5];
      qkxz = rpole[j][6];
      qkyy = rpole[j][8];
      qkyz = rpole[j][9];
      qkzz = rpole[j][12];

      // intermediates involving moments and separation distance

      dir = dix*xr + diy*yr + diz*zr;
      qix = qixx*xr + qixy*yr + qixz*zr;
      qiy = qixy*xr + qiyy*yr + qiyz*zr;
      qiz = qixz*xr + qiyz*yr + qizz*zr;
      qir = qix*xr + qiy*yr + qiz*zr;
      dkr = dkx*xr + dky*yr + dkz*zr;
      qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      qky = qkxy*xr + qkyy*yr + qkyz*zr;
      qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      qkr = qkx*xr + qky*yr + qkz*zr;
      dik = dix*dkx + diy*dky + diz*dkz;
      qik = qix*qkx + qiy*qky + qiz*qkz;
      diqk = dix*qkx + diy*qky + diz*qkz;
      dkqi = dkx*qix + dky*qiy + dkz*qiz;
      qiqk = 2.0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) +
        qixx*qkxx + qiyy*qkyy + qizz*qkzz;

      // additional intermediates involving moments and distance

      dirx = diy*zr - diz*yr;
      diry = diz*xr - dix*zr;
      dirz = dix*yr - diy*xr;
      dkrx = dky*zr - dkz*yr;
      dkry = dkz*xr - dkx*zr;
      dkrz = dkx*yr - dky*xr;
      dikx = diy*dkz - diz*dky;
      diky = diz*dkx - dix*dkz;
      dikz = dix*dky - diy*dkx;
      qirx = qiz*yr - qiy*zr;
      qiry = qix*zr - qiz*xr;
      qirz = qiy*xr - qix*yr;
      qkrx = qkz*yr - qky*zr;
      qkry = qkx*zr - qkz*xr;
      qkrz = qky*xr - qkx*yr;
      qikx = qky*qiz - qkz*qiy;
      qiky = qkz*qix - qkx*qiz;
      qikz = qkx*qiy - qky*qix;
      qixk = qixx*qkx + qixy*qky + qixz*qkz;
      qiyk = qixy*qkx + qiyy*qky + qiyz*qkz;
      qizk = qixz*qkx + qiyz*qky + qizz*qkz;
      qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz;
      qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz;
      qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz;
      qikrx = qizk*yr - qiyk*zr;
      qikry = qixk*zr - qizk*xr;
      qikrz = qiyk*xr - qixk*yr;
      qkirx = qkzi*yr - qkyi*zr;
      qkiry = qkxi*zr - qkzi*xr;
      qkirz = qkyi*xr - qkxi*yr;
      diqkx = dix*qkxx + diy*qkxy + diz*qkxz;
      diqky = dix*qkxy + diy*qkyy + diz*qkyz;
      diqkz = dix*qkxz + diy*qkyz + diz*qkzz;
      dkqix = dkx*qixx + dky*qixy + dkz*qixz;
      dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz;
      dkqiz = dkx*qixz + dky*qiyz + dkz*qizz;
      diqkrx = diqkz*yr - diqky*zr;
      diqkry = diqkx*zr - diqkz*xr;
      diqkrz = diqky*xr - diqkx*yr;
      dkqirx = dkqiz*yr - dkqiy*zr;
      dkqiry = dkqix*zr - dkqiz*xr;
      dkqirz = dkqiy*xr - dkqix*yr;
      dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy -
        2.0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz - qixz*qkxy-qiyz*qkyy-qizz*qkyz);
      dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz -
        2.0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz - qixx*qkxz-qixy*qkyz-qixz*qkzz);
      dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix -
        2.0*(qixx*qkxy+qixy*qkyy+qixz*qkyz - qixy*qkxx-qiyy*qkxy-qiyz*qkxz);

      // get reciprocal distance terms for this interaction

      rr1 = felec / r;
      rr3 = rr1 / r2;
      rr5 = 3.0 * rr3 / r2;
      rr7 = 5.0 * rr5 / r2;
      rr9 = 7.0 * rr7 / r2;
      rr11 = 9.0 * rr9 / r2;

      // calculate the real space Ewald error function terms

      ralpha = aewald * r;
      bn[0] = erfc(ralpha) / r;
      alsq2 = 2.0 * aewald*aewald;
      alsq2n = 0.0;
      if (aewald > 0.0) alsq2n = 1.0 / (MY_PIS*aewald);
      exp2a = exp(-ralpha*ralpha);
      for (k = 1; k < 6; k++) {
        bfac = (double) (k+k-1);
        alsq2n = alsq2 * alsq2n;
        bn[k] = (bfac*bn[k-1]+alsq2n*exp2a) / r2;
      }
      for (k = 0; k < 6; k++) bn[k] *= felec;

      // find damped multipole intermediates and energy value

      if (!amoeba) {
        corek = pcore[jclass];
        alphak = palpha[jclass];
        valk = pval[j];

        term1 = corei*corek;
        term1i = corek*vali;
        term2i = corek*dir;
        term3i = corek*qir;
        term1k = corei*valk;
        term2k = -corei*dkr;
        term3k = corei*qkr;
        term1ik = vali*valk;
        term2ik = valk*dir - vali*dkr + dik;
        term3ik = vali*qkr + valk*qir - dir*dkr + 2.0*(dkqi-diqk+qiqk);
        term4ik = dir*qkr - dkr*qir - 4.0*qik;
        term5ik = qir*qkr;
        damppole(r,11,alphai,alphak,dmpi,dmpj,dmpij);
        scalek = factor_mpole;
        rr1i = bn[0] - (1.0-scalek*dmpi[0])*rr1;
        rr3i = bn[1] - (1.0-scalek*dmpi[2])*rr3;
        rr5i = bn[2] - (1.0-scalek*dmpi[4])*rr5;
        rr7i = bn[3] - (1.0-scalek*dmpi[6])*rr7;
        rr1k = bn[0] - (1.0-scalek*dmpj[0])*rr1;
        rr3k = bn[1] - (1.0-scalek*dmpj[2])*rr3;
        rr5k = bn[2] - (1.0-scalek*dmpj[4])*rr5;
        rr7k = bn[3] - (1.0-scalek*dmpj[6])*rr7;
        rr1ik = bn[0] - (1.0-scalek*dmpij[0])*rr1;
        rr3ik = bn[1] - (1.0-scalek*dmpij[2])*rr3;
        rr5ik = bn[2] - (1.0-scalek*dmpij[4])*rr5;
        rr7ik = bn[3] - (1.0-scalek*dmpij[6])*rr7;
        rr9ik = bn[4] - (1.0-scalek*dmpij[8])*rr9;
        rr11ik = bn[5] - (1.0-scalek*dmpij[10])*rr11;
        rr1 = bn[0] - (1.0-scalek)*rr1;
        rr3 = bn[1] - (1.0-scalek)*rr3;
        e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik +
          term1i*rr1i + term1k*rr1k + term1ik*rr1ik +
          term2i*rr3i + term2k*rr3k + term2ik*rr3ik +
          term3i*rr5i + term3k*rr5k + term3ik*rr5ik;

        // find damped multipole intermediates for force and torque

        de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik +
          term1i*rr3i + term1k*rr3k + term1ik*rr3ik +
          term2i*rr5i + term2k*rr5k + term2ik*rr5ik +
          term3i*rr7i + term3k*rr7k + term3ik*rr7ik;
        term1 = -corek*rr3i - valk*rr3ik + dkr*rr5ik - qkr*rr7ik;
        term2 = corei*rr3k + vali*rr3ik + dir*rr5ik + qir*rr7ik;
        term3 = 2.0 * rr5ik;
        term4 = -2.0 * (corek*rr5i+valk*rr5ik - dkr*rr7ik+qkr*rr9ik);
        term5 = -2.0 * (corei*rr5k+vali*rr5ik + dir*rr7ik+qir*rr9ik);
        term6 = 4.0 * rr7ik;
        rr3 = rr3ik;

        // find standard multipole intermediates and energy value

      } else {
        term1 = ci*ck;
        term2 = ck*dir - ci*dkr + dik;
        term3 = ci*qkr + ck*qir - dir*dkr + 2.0*(dkqi-diqk+qiqk);
        term4 = dir*qkr - dkr*qir - 4.0*qik;
        term5 = qir*qkr;
        scalek = 1.0 - factor_mpole;
        rr1 = bn[0] - scalek*rr1;
        rr3 = bn[1] - scalek*rr3;
        rr5 = bn[2] - scalek*rr5;
        rr7 = bn[3] - scalek*rr7;
        rr9 = bn[4] - scalek*rr9;
        rr11 = bn[5] - scalek*rr11;
        e = term1*rr1 + term2*rr3 + term3*rr5 + term4*rr7 + term5*rr9;

        // find standard multipole intermediates for force and torque

        de = term1*rr3 + term2*rr5 + term3*rr7 + term4*rr9 + term5*rr11;
        term1 = -ck*rr3 + dkr*rr5 - qkr*rr7;
        term2 = ci*rr3 + dir*rr5 + qir*rr7;
        term3 = 2.0 * rr5;
        term4 = 2.0 * (-ck*rr5+dkr*rr7-qkr*rr9);
        term5 = 2.0 * (-ci*rr5-dir*rr7-qir*rr9);
        term6 = 4.0 * rr7;
      }

      empole += e;

      // compute the force components for this interaction

      frcx = de*xr + term1*dix + term2*dkx + term3*(diqkx-dkqix) +
        term4*qix + term5*qkx + term6*(qixk+qkxi);
      frcy = de*yr + term1*diy + term2*dky + term3*(diqky-dkqiy) +
        term4*qiy + term5*qky + term6*(qiyk+qkyi);
      frcz = de*zr + term1*diz + term2*dkz + term3*(diqkz-dkqiz) +
        term4*qiz + term5*qkz + term6*(qizk+qkzi);

      // compute the torque components for this interaction

      ttmi[0] = -rr3*dikx + term1*dirx + term3*(dqikx+dkqirx) -
        term4*qirx - term6*(qikrx+qikx);
      ttmi[1] = -rr3*diky + term1*diry + term3*(dqiky+dkqiry) -
        term4*qiry - term6*(qikry+qiky);
      ttmi[2] = -rr3*dikz + term1*dirz + term3*(dqikz+dkqirz) -
        term4*qirz - term6*(qikrz+qikz);
      ttmk[0] = rr3*dikx + term2*dkrx - term3*(dqikx+diqkrx) -
        term5*qkrx - term6*(qkirx-qikx);
      ttmk[1] = rr3*diky + term2*dkry - term3*(dqiky+diqkry) -
        term5*qkry - term6*(qkiry-qiky);
      ttmk[2] = rr3*dikz + term2*dkrz - term3*(dqikz+diqkrz) -
        term5*qkrz - term6*(qkirz-qikz);

      // increment force-based gradient and torque on first site

      f[i][0] -= frcx;
      f[i][1] -= frcy;
      f[i][2] -= frcz;
      tq[i][0] += ttmi[0];
      tq[i][1] += ttmi[1];
      tq[i][2] += ttmi[2];

      // increment force-based gradient and torque on second site

      f[j][0] += frcx;
      f[j][1] += frcy;
      f[j][2] += frcz;
      tq[j][0] += ttmk[0];
      tq[j][1] += ttmk[1];
      tq[j][2] += ttmk[2];

      // increment the virial due to pairwise Cartesian forces

      if (vflag_global) {
        vxx = -xr * frcx;
        vxy = -0.5 * (yr*frcx+xr*frcy);
        vxz = -0.5 * (zr*frcx+xr*frcz);
        vyy = -yr * frcy;
        vyz = -0.5 * (zr*frcy+yr*frcz);
        vzz = -zr * frcz;

        virmpole[0] -= vxx;
        virmpole[1] -= vyy;
        virmpole[2] -= vzz;
        virmpole[3] -= vxy;
        virmpole[4] -= vxz;
        virmpole[5] -= vyz;
      }
    }
  }

  // reverse comm to sum torque from ghost atoms to owned atoms

  crstyle = TORQUE;
  comm->reverse_comm(this);

  // resolve site torques then increment forces and virial

  for (i = 0; i < nlocal; i++) {
    torque2force(i,tq[i],fix,fiy,fiz,f);

    if (!vflag_global) continue;

    iz = zaxis2local[i];
    ix = xaxis2local[i];
    iy = yaxis2local[i];

    xiz = x[iz][0] - x[i][0];
    yiz = x[iz][1] - x[i][1];
    ziz = x[iz][2] - x[i][2];
    xix = x[ix][0] - x[i][0];
    yix = x[ix][1] - x[i][1];
    zix = x[ix][2] - x[i][2];
    xiy = x[iy][0] - x[i][0];
    yiy = x[iy][1] - x[i][1];
    ziy = x[iy][2] - x[i][2];

    vxx = xix*fix[0] + xiy*fiy[0] + xiz*fiz[0];
    vxy = 0.5 * (yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] +
                 xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
    vxz = 0.5 * (zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] +
                 xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
    vyy = yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
    vyz = 0.5 * (zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] +
                 yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);
    vzz = zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];

    virmpole[0] -= vxx;
    virmpole[1] -= vyy;
    virmpole[2] -= vzz;
    virmpole[3] -= vxy;
    virmpole[4] -= vxz;
    virmpole[5] -= vyz;
  }
}

/* ----------------------------------------------------------------------
   multipole_kspace = KSpace portion of multipole interactions
   adapted from Tinker emrecip1() routine
   literature reference:
   C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
   Representation of Electrostatics in Classical Force Fields:
   Efficient Implementation of Multipolar Interactions in
   Biomolecular Simulations", Journal of Chemical Physics, 120,
   73-87 (2004)
------------------------------------------------------------------------- */

void PairAmoeba::multipole_kspace()
{
  int i,j,k,n,ix,iy,iz;
  int nhalf1,nhalf2,nhalf3;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  double e,eterm,felec;
  double r1,r2,r3;
  double h1,h2,h3;
  double f1,f2,f3;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz,vxy,vxz,vyz;
  double volterm,denom;
  double hsq,expterm;
  double term,pterm;
  double vterm,struc2;
  double tem[3],fix[3],fiy[3],fiz[3];

  // indices into the electrostatic field array
  // decremented by 1 versus Fortran

  int deriv1[10] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
  int deriv2[10] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
  int deriv3[10] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};

  // return if the Ewald coefficient is zero

  if (aewald < 1.0e-6) return;

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  double volbox = domain->prd[0] * domain->prd[1] * domain->prd[2];

  felec = electric / am_dielectric;

  // FFT moduli pre-computations
  // set igrid for each atom and its B-spline coeffs

  nfft1 = m_kspace->nx;
  nfft2 = m_kspace->ny;
  nfft3 = m_kspace->nz;
  bsorder = m_kspace->order;

  moduli();
  bspline_fill();

  // copy multipole info to Cartesian cmp

  for (i = 0; i < nlocal; i++) {
    cmp[i][0] = rpole[i][0];
    cmp[i][1] = rpole[i][1];
    cmp[i][2] = rpole[i][2];
    cmp[i][3] = rpole[i][3];
    cmp[i][4] = rpole[i][4];
    cmp[i][5] = rpole[i][8];
    cmp[i][6] = rpole[i][12];
    cmp[i][7] = 2.0 * rpole[i][5];
    cmp[i][8] = 2.0 * rpole[i][6];
    cmp[i][9] = 2.0 * rpole[i][9];
  }

  // convert Cartesian multipoles to fractional multipoles

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values

  double ***gridpre = (double ***) m_kspace->zero();

  // map atoms to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

  double *gridfft = m_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  // zero virial accumulation variables

  vxx = vyy = vzz = vxy = vxz = vyz = 0.0;

  // perform convolution on K-space points I own

  nhalf1 = (nfft1+1) / 2;
  nhalf2 = (nfft2+1) / 2;
  nhalf3 = (nfft3+1) / 2;

  nxlo = m_kspace->nxlo_fft;
  nxhi = m_kspace->nxhi_fft;
  nylo = m_kspace->nylo_fft;
  nyhi = m_kspace->nyhi_fft;
  nzlo = m_kspace->nzlo_fft;
  nzhi = m_kspace->nzhi_fft;

  pterm = square(MY_PI/aewald);
  volterm = MY_PI * volbox;

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
        term = -pterm * hsq;
        expterm = 0.0;
        if (term > -50.0 && hsq != 0.0) {
          denom = volterm*hsq*bsmod1[i]*bsmod2[j]*bsmod3[k];
          expterm = exp(term) / denom;
          struc2 = gridfft[n]*gridfft[n] + gridfft[n+1]*gridfft[n+1];
          eterm = 0.5 * felec * expterm * struc2;
          vterm = (2.0/hsq) * (1.0-term) * eterm;
          vxx += h1*h1*vterm - eterm;
          vyy += h2*h2*vterm - eterm;
          vzz += h3*h3*vterm - eterm;
          vxy += h1*h2*vterm;
          vxz += h1*h3*vterm;
          vyz += h2*h3*vterm;
        }
        gridfft[n] *= expterm;
        gridfft[n+1] *= expterm;
        n += 2;
      }
    }
  }

  // save multipole virial for use in polarization computation

  vmsave[0] = vxx;
  vmsave[1] = vyy;
  vmsave[2] = vzz;
  vmsave[3] = vxy;
  vmsave[4] = vxz;
  vmsave[5] = vyz;

  // post-convolution operations including backward FFT
  // gridppost = my portion of 3d grid in brick decomp w/ ghost values

  double ***gridpost = (double ***) m_kspace->post_convolution();

  // get potential

  fphi_mpole(gridpost,fphi);

  for (i = 0; i < nlocal; i++) {
    for (k = 0; k < 20; k++)
      fphi[i][k] *= felec;
  }

  // convert field from fractional to Cartesian

  fphi_to_cphi(fphi,cphi);

  // increment the permanent multipole energy and gradient

  e = 0.0;
  for (i = 0; i < nlocal; i++) {
    f1 = 0.0;
    f2 = 0.0;
    f3 = 0.0;
    for (k = 0; k < 10; k++) {
      e += fmp[i][k]*fphi[i][k];
      f1 += fmp[i][k]*fphi[i][deriv1[k]];
      f2 += fmp[i][k]*fphi[i][deriv2[k]];
      f3 += fmp[i][k]*fphi[i][deriv3[k]];
    }
    f1 *= nfft1;
    f2 *= nfft2;
    f3 *= nfft3;
    h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;  // matvec?
    h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
    h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;
    f[i][0] -= h1;
    f[i][1] -= h2;
    f[i][2] -= h3;
  }
  empole += 0.5*e;

  // augment the permanent multipole virial contributions

  if (vflag_global) {
    for (i = 0; i < nlocal; i++) {
      vxx = vxx - cmp[i][1]*cphi[i][1] - 2.0*cmp[i][4]*cphi[i][4] -
        cmp[i][7]*cphi[i][7] - cmp[i][8]*cphi[i][8];
      vxy = vxy - 0.5*(cmp[i][2]*cphi[i][1]+cmp[i][1]*cphi[i][2]) -
        (cmp[i][4]+cmp[i][5])*cphi[i][7] - 0.5*cmp[i][7]*(cphi[i][4]+cphi[i][5]) -
        0.5*(cmp[i][8]*cphi[i][9]+cmp[i][9]*cphi[i][8]);
      vxz = vxz - 0.5*(cmp[i][3]*cphi[i][1]+cmp[i][1]*cphi[i][3]) -
        (cmp[i][4]+cmp[i][6])*cphi[i][8] - 0.5*cmp[i][8]*(cphi[i][4]+cphi[i][6]) -
        0.5*(cmp[i][7]*cphi[i][9]+cmp[i][9]*cphi[i][7]);
      vyy = vyy - cmp[i][2]*cphi[i][2] - 2.0*cmp[i][5]*cphi[i][5] -
        cmp[i][7]*cphi[i][7] - cmp[i][9]*cphi[i][9];
      vyz = vyz - 0.5*(cmp[i][3]*cphi[i][2]+cmp[i][2]*cphi[i][3]) -
        (cmp[i][5]+cmp[i][6])*cphi[i][9] - 0.5*cmp[i][9]*(cphi[i][5]+cphi[i][6]) -
        0.5*(cmp[i][7]*cphi[i][8]+cmp[i][8]*cphi[i][7]);
      vzz = vzz - cmp[i][3]*cphi[i][3] - 2.0*cmp[i][6]*cphi[i][6] -
        cmp[i][8]*cphi[i][8] - cmp[i][9]*cphi[i][9];
    }
  }

  // resolve site torques then increment forces and virial

  for (i = 0; i < nlocal; i++) {
    tem[0] = cmp[i][3]*cphi[i][2] - cmp[i][2]*cphi[i][3] +
      2.0*(cmp[i][6]-cmp[i][5])*cphi[i][9] +
      cmp[i][8]*cphi[i][7] + cmp[i][9]*cphi[i][5] -
      cmp[i][7]*cphi[i][8] - cmp[i][9]*cphi[i][6];
    tem[1] = cmp[i][1]*cphi[i][3] - cmp[i][3]*cphi[i][1] +
      2.0*(cmp[i][4]-cmp[i][6])*cphi[i][8] +
      cmp[i][7]*cphi[i][9] + cmp[i][8]*cphi[i][6] -
        cmp[i][8]*cphi[i][4] - cmp[i][9]*cphi[i][7];
    tem[2] = cmp[i][2]*cphi[i][1] - cmp[i][1]*cphi[i][2] +
      2.0*(cmp[i][5]-cmp[i][4])*cphi[i][7] +
      cmp[i][7]*cphi[i][4] + cmp[i][9]*cphi[i][8] -
      cmp[i][7]*cphi[i][5] - cmp[i][8]*cphi[i][9];

    torque2force(i,tem,fix,fiy,fiz,f);

    if (vflag_global) {
      iz = zaxis2local[i];
      ix = xaxis2local[i];
      iy = yaxis2local[i];

      xiz = x[iz][0] - x[i][0];
      yiz = x[iz][1] - x[i][1];
      ziz = x[iz][2] - x[i][2];
      xix = x[ix][0] - x[i][0];
      yix = x[ix][1] - x[i][1];
      zix = x[ix][2] - x[i][2];
      xiy = x[iy][0] - x[i][0];
      yiy = x[iy][1] - x[i][1];
      ziy = x[iy][2] - x[i][2];

      vxx += xix*fix[0] + xiy*fiy[0] + xiz*fiz[0];
      vxy += 0.5*(yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] +
                  xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
      vxz += 0.5*(zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] +
                  xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
      vyy += yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
      vyz += 0.5*(zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] +
                  yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);
      vzz += zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];
    }
  }

  // increment total internal virial tensor components

  if (vflag_global) {
    virmpole[0] -= vxx;
    virmpole[1] -= vyy;
    virmpole[2] -= vzz;
    virmpole[3] -= vxy;
    virmpole[4] -= vxz;
    virmpole[5] -= vyz;
  }
}

/* ----------------------------------------------------------------------
   damppole generates coefficients for the charge penetration
   damping function for powers of the interatomic distance

   literature references:

   L. V. Slipchenko and M. S. Gordon, "Electrostatic Energy in the
   Effective Fragment Potential Method: Theory and Application to
   the Benzene Dimer", Journal of Computational Chemistry, 28,
   276-291 (2007)  [Gordon f1 and f2 models]

   J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and
   J. W. Ponder, "An Optimized Charge Penetration Model for Use with
   the AMOEBA Force Field", Physical Chemistry Chemical Physics, 19,
   276-291 (2017)
------------------------------------------------------------------------- */

void PairAmoeba::damppole(double r, int rorder, double alphai, double alphak,
                          double *dmpi, double *dmpk, double *dmpik)
{
  double termi,termk;
  double termi2,termk2;
  double alphai2,alphak2;
  double eps,diff;
  double expi,expk;
  double dampi,dampk;
  double dampi2,dampi3;
  double dampi4,dampi5;
  double dampi6,dampi7;
  double dampi8;
  double dampk2,dampk3;
  double dampk4,dampk5;
  double dampk6;

  // compute tolerance and exponential damping factors

  eps = 0.001;
  diff = fabs(alphai-alphak);
  dampi = alphai * r;
  dampk = alphak * r;
  expi = exp(-dampi);
  expk = exp(-dampk);

  // core-valence charge penetration damping for Gordon f1

  dampi2 = dampi * dampi;
  dampi3 = dampi * dampi2;
  dampi4 = dampi2 * dampi2;
  dampi5 = dampi2 * dampi3;
  dmpi[0] = 1.0 - (1.0 + 0.5*dampi)*expi;
  dmpi[2] = 1.0 - (1.0 + dampi + 0.5*dampi2)*expi;
  dmpi[4] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0)*expi;
  dmpi[6] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 + dampi4/30.0)*expi;
  dmpi[8] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                   4.0*dampi4/105.0 + dampi5/210.0)*expi;
  if (diff < eps) {
    dmpk[0] = dmpi[0];
    dmpk[2] = dmpi[2];
    dmpk[4] = dmpi[4];
    dmpk[6] = dmpi[6];
    dmpk[8] = dmpi[8];
  } else {
    dampk2 = dampk * dampk;
    dampk3 = dampk * dampk2;
    dampk4 = dampk2 * dampk2;
    dampk5 = dampk2 * dampk3;
    dmpk[0] = 1.0 - (1.0 + 0.5*dampk)*expk;
    dmpk[2] = 1.0 - (1.0 + dampk + 0.5*dampk2)*expk;
    dmpk[4] = 1.0 - (1.0 + dampk + 0.5*dampk2 + dampk3/6.0)*expk;
    dmpk[6] = 1.0 - (1.0 + dampk + 0.5*dampk2 + dampk3/6.0 + dampk4/30.0)*expk;
    dmpk[8] = 1.0 - (1.0 + dampk + 0.5*dampk2 + dampk3/6.0 +
                     4.0*dampk4/105.0 + dampk5/210.0)*expk;
  }

  // valence-valence charge penetration damping for Gordon f1

  if (diff < eps) {
    dampi6 = dampi3 * dampi3;
    dampi7 = dampi3 * dampi4;
    dmpik[0] = 1.0 - (1.0 + 11.0*dampi/16.0 + 3.0*dampi2/16.0 +
                      dampi3/48.0)*expi;
    dmpik[2] = 1.0 - (1.0 + dampi + 0.5*dampi2 +
                      7.0*dampi3/48.0 + dampi4/48.0)*expi;
    dmpik[4] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                      dampi4/24.0 + dampi5/144.0)*expi;
    dmpik[6] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                      dampi4/24.0 + dampi5/120.0 + dampi6/720.0)*expi;
    dmpik[8] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                      dampi4/24.0 + dampi5/120.0 + dampi6/720.0 +
                      dampi7/5040.0)*expi;
    if (rorder >= 11) {
      dampi8 = dampi4 * dampi4;
      dmpik[10] = 1.0 - (1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                         dampi4/24.0 + dampi5/120.0 + dampi6/720.0 +
                         dampi7/5040.0 + dampi8/45360.0)*expi;
    }

  } else {
    alphai2 = alphai * alphai;
    alphak2 = alphak * alphak;
    termi = alphak2 / (alphak2-alphai2);
    termk = alphai2 / (alphai2-alphak2);
    termi2 = termi * termi;
    termk2 = termk * termk;
    dmpik[0] = 1.0 - termi2*(1.0 + 2.0*termk + 0.5*dampi)*expi -
      termk2*(1.0 + 2.0*termi + 0.5*dampk)*expk;
    dmpik[2] = 1.0 - termi2*(1.0+dampi+0.5*dampi2)*expi -
      termk2*(1.0+dampk+0.5*dampk2)*expk -
      2.0*termi2*termk*(1.0+dampi)*expi -
      2.0*termk2*termi*(1.0+dampk)*expk;
    dmpik[4] = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2 + dampi3/6.0)*expi -
      termk2*(1.0 + dampk + 0.5*dampk2 + dampk3/6.0)*expk -
      2.0*termi2*termk*(1.0 + dampi + dampi2/3.0)*expi -
      2.0*termk2*termi*(1.0 + dampk + dampk2/3.0)*expk;
    dmpik[6] = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2 +
                             dampi3/6.0 + dampi4/30.0)*expi -
      termk2*(1.0 + dampk + 0.5*dampk2 + dampk3/6.0 + dampk4/30.0)*expk -
      2.0*termi2*termk*(1.0 + dampi + 2.0*dampi2/5.0 + dampi3/15.0)*expi -
      2.0*termk2*termi*(1.0 + dampk + 2.0*dampk2/5.0 + dampk3/15.0)*expk;
    dmpik[8] = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                             4.0*dampi4/105.0 + dampi5/210.0)*expi -
      termk2*(1.0 + dampk + 0.5*dampk2 + dampk3/6.0 +
              4.0*dampk4/105.0 + dampk5/210.0)*expk -
      2.0*termi2*termk*(1.0 + dampi + 3.0*dampi2/7.0 +
                        2.0*dampi3/21.0 + dampi4/105.0)*expi -
      2.0*termk2*termi*(1.0 + dampk + 3.0*dampk2/7.0 +
                        2.0*dampk3/21.0 + dampk4/105.0)*expk;

    if (rorder >= 11) {
      dampi6 = dampi3 * dampi3;
      dampk6 = dampk3 * dampk3;
      dmpik[10] = 1.0 - termi2*(1.0 + dampi + 0.5*dampi2 + dampi3/6.0 +
                                5.0*dampi4/126.0 + 2.0*dampi5/315.0 +
                                dampi6/1890.0)*expi -
        termk2*(1.0 + dampk + 0.5*dampk2 + dampk3/6.0 + 5.0*dampk4/126.0 +
                2.0*dampk5/315.0 + dampk6/1890.0)*expk -
        2.0*termi2*termk*(1.0 + dampi + 4.0*dampi2/9.0 + dampi3/9.0 +
                          dampi4/63.0 + dampi5/945.0)*expi -
        2.0*termk2*termi*(1.0 + dampk + 4.0*dampk2/9.0 + dampk3/9.0 +
                          dampk4/63.0 + dampk5/945.0)*expk;
    }
  }
}
