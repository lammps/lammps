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

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>

using namespace LAMMPS_NS;

enum{FIELD,ZRSD,TORQUE,UFLD};                          // reverse comm
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   repulsion = Pauli repulsion interactions
   adapted from Tinker erepel1b() routine
------------------------------------------------------------------------- */

void PairAmoeba::repulsion()
{
  int i,j,k,ii,jj,itype,jtype;
  int ix,iy,iz;
  double e;
  double eterm,de;
  double xi,yi,zi;
  double xr,yr,zr;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double r,r2,r3,r4,r5;
  double rr1,rr3,rr5;
  double rr7,rr9,rr11;
  double dix,diy,diz;
  double qixx,qixy,qixz;
  double qiyy,qiyz,qizz;
  double dkx,dky,dkz;
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
  double term1,term2,term3;
  double term4,term5,term6;
  double sizi,sizk,sizik;
  double vali,valk;
  double dmpi,dmpk;
  double frcx,frcy,frcz;
  double taper,dtaper;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double factor_repel;
  double ttri[3],ttrk[3];
  double fix[3],fiy[3],fiz[3];
  double dmpik[11];

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // set cutoffs and taper coeffs

  choose(REPULSE);

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  // zero repulsion torque on owned + ghost atoms

  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    tq[i][0] = 0.0;
    tq[i][1] = 0.0;
    tq[i][2] = 0.0;
  }

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // double loop over owned atoms and neighbors

  // DEBUG
  //FILE *fp = fopen("lammps.dat","w");

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];
    sizi = sizpr[itype];
    dmpi = dmppr[itype];
    vali = elepr[itype];
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_repel = special_repel[sbmask15(j)];
      if (factor_repel == 0.0) continue;
      j &= NEIGHMASK15;

      xr = x[j][0] - xi;
      yr = x[j][1] - yi;
      zr = x[j][2] - zi;
      r2 = xr*xr + yr*yr + zr*zr;
      if (r2 > off2) continue;

      jtype = amtype[j];

      r = sqrt(r2);
      sizk = sizpr[jtype];
      dmpk = dmppr[jtype];
      valk = elepr[jtype];
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
        2.0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz-qixz*qkxy-qiyz*qkyy-qizz*qkyz);
      dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz -
        2.0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz-qixx*qkxz-qixy*qkyz-qixz*qkzz);
      dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix -
        2.0*(qixx*qkxy+qixy*qkyy+qixz*qkyz-qixy*qkxx-qiyy*qkxy-qiyz*qkxz);

      // get reciprocal distance terms for this interaction

      rr1 = 1.0 / r;
      rr3 = rr1 / r2;
      rr5 = 3.0 * rr3 / r2;
      rr7 = 5.0 * rr5 / r2;
      rr9 = 7.0 * rr7 / r2;
      rr11 = 9.0 * rr9 / r2;

      // get damping coefficients for the Pauli repulsion energy

      damprep(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,11,dmpi,dmpk,dmpik);

      // calculate intermediate terms needed for the energy

      term1 = vali*valk;
      term2 = valk*dir - vali*dkr + dik;
      term3 = vali*qkr + valk*qir - dir*dkr + 2.0*(dkqi-diqk+qiqk);
      term4 = dir*qkr - dkr*qir - 4.0*qik;
      term5 = qir*qkr;
      eterm = term1*dmpik[0] + term2*dmpik[2] +
        term3*dmpik[4] + term4*dmpik[6] + term5*dmpik[8];

      // compute the Pauli repulsion energy for this interaction

      sizik = sizi * sizk * factor_repel;
      e = sizik * eterm * rr1;

      // calculate intermediate terms for force and torque

      de = term1*dmpik[2] + term2*dmpik[4] + term3*dmpik[6] +
        term4*dmpik[8] + term5*dmpik[10];
      term1 = -valk*dmpik[2] + dkr*dmpik[4] - qkr*dmpik[6];
      term2 = vali*dmpik[2] + dir*dmpik[4] + qir*dmpik[6];
      term3 = 2.0 * dmpik[4];
      term4 = 2.0 * (-valk*dmpik[4] + dkr*dmpik[6] - qkr*dmpik[8]);
      term5 = 2.0 * (-vali*dmpik[4] - dir*dmpik[6] - qir*dmpik[8]);
      term6 = 4.0 * dmpik[6];

      // compute the force components for this interaction

      frcx = de*xr + term1*dix + term2*dkx + term3*(diqkx-dkqix) +
        term4*qix + term5*qkx + term6*(qixk+qkxi);
      frcy = de*yr + term1*diy + term2*dky + term3*(diqky-dkqiy) +
        term4*qiy + term5*qky + term6*(qiyk+qkyi);
      frcz = de*zr + term1*diz + term2*dkz + term3*(diqkz-dkqiz) +
        term4*qiz + term5*qkz + term6*(qizk+qkzi);
      frcx = frcx*rr1 + eterm*rr3*xr;
      frcy = frcy*rr1 + eterm*rr3*yr;
      frcz = frcz*rr1 + eterm*rr3*zr;
      frcx = sizik * frcx;
      frcy = sizik * frcy;
      frcz = sizik * frcz;

      // compute the torque components for this interaction

      ttri[0] = -dmpik[2]*dikx + term1*dirx + term3*(dqikx+dkqirx) -
        term4*qirx - term6*(qikrx+qikx);
      ttri[1] = -dmpik[2]*diky + term1*diry + term3*(dqiky+dkqiry) -
        term4*qiry - term6*(qikry+qiky);
      ttri[2] = -dmpik[2]*dikz + term1*dirz + term3*(dqikz+dkqirz) -
        term4*qirz - term6*(qikrz+qikz);
      ttrk[0] = dmpik[2]*dikx + term2*dkrx - term3*(dqikx+diqkrx) -
        term5*qkrx - term6*(qkirx-qikx);
      ttrk[1] = dmpik[2]*diky + term2*dkry - term3*(dqiky+diqkry) -
        term5*qkry - term6*(qkiry-qiky);
      ttrk[2] = dmpik[2]*dikz + term2*dkrz - term3*(dqikz+diqkrz) -
        term5*qkrz - term6*(qkirz-qikz);
      ttri[0] = sizik * ttri[0] * rr1;
      ttri[1] = sizik * ttri[1] * rr1;
      ttri[2] = sizik * ttri[2] * rr1;
      ttrk[0] = sizik * ttrk[0] * rr1;
      ttrk[1] = sizik * ttrk[1] * rr1;
      ttrk[2] = sizik * ttrk[2] * rr1;

      // use energy switching if near the cutoff distance

      if (r2 > cut2) {
        r3 = r2 * r;
        r4 = r2 * r2;
        r5 = r2 * r3;
        taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0;
        dtaper = 5.0*c5*r4 + 4.0*c4*r3 + 3.0*c3*r2 + 2.0*c2*r + c1;
        dtaper *= e * rr1;
        e *= taper;
        frcx = frcx*taper - dtaper*xr;
        frcy = frcy*taper - dtaper*yr;
        frcz = frcz*taper - dtaper*zr;
        for (k = 0; k < 3; k++) {
          ttri[k] *= taper;
          ttrk[k] *= taper;
        }
      }

      erepulse += e;

      // increment force-based gradient and torque on atom I

      f[i][0] -= frcx;
      f[i][1] -= frcy;
      f[i][2] -= frcz;
      tq[i][0] += ttri[0];
      tq[i][1] += ttri[1];
      tq[i][2] += ttri[2];

      // increment force-based gradient and torque on atom J

      f[j][0] += frcx;
      f[j][1] += frcy;
      f[j][2] += frcz;
      tq[j][0] += ttrk[0];
      tq[j][1] += ttrk[1];
      tq[j][2] += ttrk[2];

      // increment the virial due to pairwise Cartesian forces

      if (vflag_global) {
        vxx = -xr * frcx;
        vxy = -0.5 * (yr*frcx+xr*frcy);
        vxz = -0.5 * (zr*frcx+xr*frcz);
        vyy = -yr * frcy;
        vyz = -0.5 * (zr*frcy+yr*frcz);
        vzz = -zr * frcz;

        virrepulse[0] -= vxx;
        virrepulse[1] -= vyy;
        virrepulse[2] -= vzz;
        virrepulse[3] -= vxy;
        virrepulse[4] -= vxz;
        virrepulse[5] -= vyz;
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
    vyy = yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
    vzz = zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];
    vxy = 0.5 * (yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] +
                 xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
    vxz = 0.5 * (zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] +
                 xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
    vyz = 0.5 * (zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] +
                 yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);

    virrepulse[0] -= vxx;
    virrepulse[1] -= vyy;
    virrepulse[2] -= vzz;
    virrepulse[3] -= vxy;
    virrepulse[4] -= vxz;
    virrepulse[5] -= vyz;
  }
}

/* ----------------------------------------------------------------------
   damprep generates coefficients for the Pauli repulsion
   damping function for powers of the interatomic distance

   literature reference:

   J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion: An
   Anisotropic, Atomic Multipole Model", Journal of Chemical Physics,
   150, 084104 (2019)
------------------------------------------------------------------------- */

void PairAmoeba::damprep(double r, double r2, double rr1, double rr3,
                         double rr5, double rr7, double rr9, double rr11,
                         int rorder, double dmpi, double dmpk, double *dmpik)
{
  double r3,r4;
  double r5,r6,r7,r8;
  double s,ds,d2s;
  double d3s,d4s,d5s;
  double dmpi2,dmpk2;
  double dmpi22,dmpi23;
  double dmpi24,dmpi25;
  double dmpi26,dmpi27;
  double dmpk22,dmpk23;
  double dmpk24,dmpk25;
  double dmpk26;
  double eps,diff;
  double expi,expk;
  double dampi,dampk;
  double pre,term,tmp;

  // compute tolerance value for damping exponents

  eps = 0.001;
  diff = fabs(dmpi-dmpk);

  // treat the case where alpha damping exponents are equal

  if (diff < eps) {
    r3 = r2 * r;
    r4 = r3 * r;
    r5 = r4 * r;
    r6 = r5 * r;
    r7 = r6 * r;
    dmpi2 = 0.5 * dmpi;
    dampi = dmpi2 * r;
    expi = exp(-dampi);
    dmpi22 = dmpi2 * dmpi2;
    dmpi23 = dmpi22 * dmpi2;
    dmpi24 = dmpi23 * dmpi2;
    dmpi25 = dmpi24 * dmpi2;
    dmpi26 = dmpi25 * dmpi2;
    pre = 128.0;
    s = (r + dmpi2*r2 + dmpi22*r3/3.0) * expi;

    ds = (dmpi22*r3 + dmpi23*r4) * expi / 3.0;
    d2s = dmpi24 * expi * r5 / 9.0;
    d3s = dmpi25 * expi * r6 / 45.0;
    d4s = (dmpi25*r6 + dmpi26*r7) * expi / 315.0;
    if (rorder >= 11) {
      r8 = r7 * r;
      dmpi27 = dmpi2 * dmpi26;
      d5s = (dmpi25*r6 + dmpi26*r7 + dmpi27*r8/3.0) * expi / 945.0;
    } else d5s = 0.0;

  // treat the case where alpha damping exponents are unequal

  } else {
    r3 = r2 * r;
    r4 = r3 * r;
    r5 = r4 * r;
    dmpi2 = 0.5 * dmpi;
    dmpk2 = 0.5 * dmpk;
    dampi = dmpi2 * r;
    dampk = dmpk2 * r;
    expi = exp(-dampi);
    expk = exp(-dampk);
    dmpi22 = dmpi2 * dmpi2;
    dmpi23 = dmpi22 * dmpi2;
    dmpi24 = dmpi23 * dmpi2;
    dmpi25 = dmpi24 * dmpi2;
    dmpk22 = dmpk2 * dmpk2;
    dmpk23 = dmpk22 * dmpk2;
    dmpk24 = dmpk23 * dmpk2;
    dmpk25 = dmpk24 * dmpk2;
    term = dmpi22 - dmpk22;
    pre = 8192.0 * dmpi23 * dmpk23 / (term*term*term*term);
    tmp = 4.0 * dmpi2 * dmpk2 / term;
    s = (dampi-tmp)*expk + (dampk+tmp)*expi;

    ds = (dmpi2*dmpk2*r2 - 4.0*dmpi2*dmpk22*r/term -
          4.0*dmpi2*dmpk2/term) * expk +
      (dmpi2*dmpk2*r2 + 4.0*dmpi22*dmpk2*r/term + 4.0*dmpi2*dmpk2/term) * expi;
    d2s = (dmpi2*dmpk2*r2/3.0 + dmpi2*dmpk22*r3/3.0 -
           (4.0/3.0)*dmpi2*dmpk23*r2/term - 4.0*dmpi2*dmpk22*r/term -
           4.0*dmpi2*dmpk2/term) * expk +
      (dmpi2*dmpk2*r2/3.0 + dmpi22*dmpk2*r3/3.0 +
       (4.0/3.0)*dmpi23*dmpk2*r2/term + 4.0*dmpi22*dmpk2*r/term +
       4.0*dmpi2*dmpk2/term) * expi;
    d3s = (dmpi2*dmpk23*r4/15.0 + dmpi2*dmpk22*r3/5.0 + dmpi2*dmpk2*r2/5.0 -
           (4.0/15.0)*dmpi2*dmpk24*r3/term - (8.0/5.0)*dmpi2*dmpk23*r2/term -
           4.0*dmpi2*dmpk22*r/term - 4.0/term*dmpi2*dmpk2) * expk +
      (dmpi23*dmpk2*r4/15.0 + dmpi22*dmpk2*r3/5.0 + dmpi2*dmpk2*r2/5.0 +
       (4.0/15.0)*dmpi24*dmpk2*r3/term + (8.0/5.0)*dmpi23*dmpk2*r2/term +
       4.0*dmpi22*dmpk2*r/term + 4.0/term*dmpi2*dmpk2) * expi;
    d4s = (dmpi2*dmpk24*r5/105.0 + (2.0/35.0)*dmpi2*dmpk23*r4 +
           dmpi2*dmpk22*r3/7.0 + dmpi2*dmpk2*r2/7.0 -
           (4.0/105.0)*dmpi2*dmpk25*r4/term - (8.0/21.0)*dmpi2*dmpk24*r3/term -
           (12.0/7.0)*dmpi2*dmpk23*r2/term - 4.0*dmpi2*dmpk22*r/term -
           4.0*dmpi2*dmpk2/term) * expk +
      (dmpi24*dmpk2*r5/105.0 + (2.0/35.0)*dmpi23*dmpk2*r4 +
       dmpi22*dmpk2*r3/7.0 + dmpi2*dmpk2*r2/7.0 +
       (4.0/105.0)*dmpi25*dmpk2*r4/term + (8.0/21.0)*dmpi24*dmpk2*r3/term +
       (12.0/7.0)*dmpi23*dmpk2*r2/term + 4.0*dmpi22*dmpk2*r/term +
       4.0*dmpi2*dmpk2/term) * expi;

    if (rorder >= 11) {
      r6 = r5 * r;
      dmpi26 = dmpi25 * dmpi2;
      dmpk26 = dmpk25 * dmpk2;
      d5s = (dmpi2*dmpk25*r6/945.0 + (2.0/189.0)*dmpi2*dmpk24*r5 +
             dmpi2*dmpk23*r4/21.0 + dmpi2*dmpk22*r3/9.0 + dmpi2*dmpk2*r2/9.0 -
             (4.0/945.0)*dmpi2*dmpk26*r5/term -
             (4.0/63.0)*dmpi2*dmpk25*r4/term - (4.0/9.0)*dmpi2*dmpk24*r3/term -
             (16.0/9.0)*dmpi2*dmpk23*r2/term - 4.0*dmpi2*dmpk22*r/term -
             4.0*dmpi2*dmpk2/term) * expk +
        (dmpi25*dmpk2*r6/945.0 + (2.0/189.0)*dmpi24*dmpk2*r5 +
         dmpi23*dmpk2*r4/21.0 + dmpi22*dmpk2*r3/9.0 + dmpi2*dmpk2*r2/9.0 +
         (4.0/945.0)*dmpi26*dmpk2*r5/term + (4.0/63.0)*dmpi25*dmpk2*r4/term +
         (4.0/9.0)*dmpi24*dmpk2*r3/term + (16.0/9.0)*dmpi23*dmpk2*r2/term +
         4.0*dmpi22*dmpk2*r/term + 4.0*dmpi2*dmpk2/term) * expi;
    } else d5s = 0.0;
  }

  // convert partial derivatives into full derivatives

  s = s * rr1;
  ds = ds * rr3;
  d2s = d2s * rr5;
  d3s = d3s * rr7;
  d4s = d4s * rr9;
  d5s = d5s * rr11;
  dmpik[0] = 0.5 * pre * s * s;
  dmpik[2] = pre * s * ds;
  dmpik[4] = pre * (s*d2s + ds*ds);
  dmpik[6] = pre * (s*d3s + 3.0*ds*d2s);
  dmpik[8] = pre * (s*d4s + 4.0*ds*d3s + 3.0*d2s*d2s);
  if (rorder >= 11) dmpik[10] = pre * (s*d5s + 5.0*ds*d4s + 10.0*d2s*d3s);
}
