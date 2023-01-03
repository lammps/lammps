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

#include "amoeba_convolution.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "math_const.h"
#include "math_special.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

using MathSpecial::square;
using MathSpecial::cube;

enum{FIELD,ZRSD,TORQUE,UFLD};                          // reverse comm
enum{MUTUAL,OPT,TCG,DIRECT};
enum{VDWL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};

/* ----------------------------------------------------------------------
   polar = induced dipole polarization
   adapted from Tinker epolar1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::polar()
{
  int i;
  int ix,iy,iz;
  double felec,term;
  double dix,diy,diz;
  double uix,uiy,uiz;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double fix[3],fiy[3],fiz[3];
  double tep[3];

  // set cutoffs, taper coeffs, and PME params

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);

  // owned atoms

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  // set the energy unit conversion factor

  felec = electric / am_dielectric;

  // compute the total induced dipole polarization energy

  polar_energy();

  // compute the real space part of the dipole interactions

  if (polar_rspace_flag) polar_real();

  // compute the reciprocal space part of dipole interactions

  if (polar_kspace_flag) polar_kspace();

  // compute the Ewald self-energy torque and virial terms

  term = (4.0/3.0) * felec * cube(aewald) / MY_PIS;

  for (i = 0; i < nlocal; i++) {
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    uix = 0.5 * (uind[i][0]+uinp[i][0]);
    uiy = 0.5 * (uind[i][1]+uinp[i][1]);
    uiz = 0.5 * (uind[i][2]+uinp[i][2]);
    tep[0] = term * (diy*uiz-diz*uiy);
    tep[1] = term * (diz*uix-dix*uiz);
    tep[2] = term * (dix*uiy-diy*uix);

    torque2force(i,tep,fix,fiy,fiz,f);

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

    virpolar[0] -= vxx;
    virpolar[1] -= vyy;
    virpolar[2] -= vzz;
    virpolar[3] -= vxy;
    virpolar[4] -= vxz;
    virpolar[5] -= vyz;
  }
}

/* ----------------------------------------------------------------------
   polar_energy = inducded dipole polarization energy
   adapted from Tinker epolar1e() routine
------------------------------------------------------------------------- */

void PairAmoeba::polar_energy()
{
  int i,j,itype;
  double e,felec,fi;

  // owned atoms

  int nlocal = atom->nlocal;

  // set the energy unit conversion factor

  felec = -0.5 * electric / am_dielectric;

  // get polarization energy via induced dipoles times field

  for (i = 0; i < nlocal; i++) {
    itype = amtype[i];
    fi = felec / polarity[itype];
    e = 0.0;
    for (j = 0; j < 3; j++)
      e += fi*uind[i][j]*udirp[i][j];
    epolar += e;
  }
}

/* ----------------------------------------------------------------------
   polar_real = real-space portion of induced dipole polarization
   adapted from Tinker epreal1d() routine
------------------------------------------------------------------------- */

void PairAmoeba::polar_real()
{
  int i,j,k,m,ii,jj,jextra,itype,jtype,iclass,jclass,igroup,jgroup;
  int ix,iy,iz;
  double felec,bfac;
  double alsq2,alsq2n;
  double exp2a,ralpha;
  double damp,expdamp;
  double pdi,pti;
  double pgamma;
  double temp3,temp5,temp7;
  double sc3,sc5,sc7;
  double psc3,psc5,psc7;
  double dsc3,dsc5,dsc7;
  double usc3,usc5;
  double psr3,psr5,psr7;
  double dsr3,dsr5,dsr7;
  double usr5;
  double rr3core,rr5core;
  double rr3i,rr5i;
  double rr7i,rr9i;
  double rr3k,rr5k;
  double rr7k,rr9k;
  double rr5ik,rr7ik;
  double xi,yi,zi;
  double xr,yr,zr;
  double r,r2,rr1,rr3;
  double rr5,rr7,rr9;
  double ci,dix,diy,diz;
  double qixx,qixy,qixz;
  double qiyy,qiyz,qizz;
  double uix,uiy,uiz;
  double uixp,uiyp,uizp;
  double ck,dkx,dky,dkz;
  double qkxx,qkxy,qkxz;
  double qkyy,qkyz,qkzz;
  double ukx,uky,ukz;
  double ukxp,ukyp,ukzp;
  double dir,uir,uirp;
  double dkr,ukr,ukrp;
  double qix,qiy,qiz,qir;
  double qkx,qky,qkz,qkr;
  double corei,corek;
  double vali,valk;
  double alphai,alphak;
  double uirm,ukrm;
  double tuir,tukr;
  double tixx,tiyy,tizz;
  double tixy,tixz,tiyz;
  double tkxx,tkyy,tkzz;
  double tkxy,tkxz,tkyz;
  double tix3,tiy3,tiz3;
  double tix5,tiy5,tiz5;
  double tkx3,tky3,tkz3;
  double tkx5,tky5,tkz5;
  double term1,term2,term3;
  double term4,term5;
  double term6,term7;
  double term1core;
  double term1i,term2i,term3i;
  double term4i,term5i,term6i;
  double term7i,term8i;
  double term1k,term2k,term3k;
  double term4k,term5k,term6k;
  double term7k,term8k;
  double depx,depy,depz;
  double frcx,frcy,frcz;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double factor_pscale,factor_dscale,factor_uscale,factor_wscale;
  double rc3[3],rc5[3],rc7[3];
  double prc3[3],prc5[3],prc7[3];
  double drc3[3],drc5[3],drc7[3];
  double urc3[3],urc5[3],tep[3];
  double fix[3],fiy[3],fiz[3];
#if 0  // for poltyp TCG which is currently not supported
  double uax[3],uay[3],uaz[3];
  double ubx[3],uby[3],ubz[3];
  double uaxp[3],uayp[3],uazp[3];
  double ubxp[3],ubyp[3],ubzp[3];
#endif
  double dmpi[9],dmpk[9];
  double dmpik[9];
  double bn[5];

  int inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // owned atoms

  double *pval = atom->dvector[index_pval];
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // initialize ufld,dulfd to zero for owned and ghost atoms

  for (i = 0; i < nall; i++)
    for (j = 0; j < 3; j++)
      ufld[i][j] = 0.0;

  for (i = 0; i < nall; i++)
    for (j = 0; j < 6; j++)
      dufld[i][j] = 0.0;

  // neigh list

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // set the energy unit conversion factor
  // NOTE: why 1/2 ?

  felec = 0.5 * electric / am_dielectric;

  // compute the dipole polarization gradient components

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = amtype[i];
    iclass = amtype2class[itype];
    igroup = amgroup[i];
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
    uix = uind[i][0];
    uiy = uind[i][1];
    uiz = uind[i][2];
    uixp = uinp[i][0];
    uiyp = uinp[i][1];
    uizp = uinp[i][2];
#if 0  // for poltyp TCG which is currently not supported
    for (m = 0; m < tcgnab; m++) {
      uax[m] = uad[m][i][0];
      uay[m] = uad[m][i][1];
      uaz[m] = uad[m][i][2];
      uaxp[m] = uap[m][i][0];
      uayp[m] = uap[m][i][1];
      uazp[m] = uap[m][i][2];
      ubx[m] = ubd[m][i][0];
      uby[m] = ubd[m][i][1];
      ubz[m] = ubd[m][i][2];
      ubxp[m] = ubp[m][i][0];
      ubyp[m] = ubp[m][i][1];
      ubzp[m] = ubp[m][i][2];
    }
#endif
    if (amoeba) {
      pdi = pdamp[itype];
      pti = thole[itype];
    } else {
      corei = pcore[iclass];
      alphai = palpha[iclass];
      vali = pval[i];
    }

    // evaluate all sites within the cutoff distance

    for (jj = 0; jj < jnum; jj++) {
      jextra = jlist[jj];
      j = jextra & NEIGHMASK15;

      xr = x[j][0] - xi;
      yr = x[j][1] - yi;
      zr = x[j][2] - zi;
      r2 = xr*xr + yr*yr + zr*zr;
      if (r2 > off2) continue;

      jtype = amtype[j];
      jclass = amtype2class[jtype];
      jgroup = amgroup[j];

      if (amoeba) {
        factor_wscale = special_polar_wscale[sbmask15(jextra)];
        if (igroup == jgroup) {
          factor_pscale = special_polar_piscale[sbmask15(jextra)];
          factor_dscale = polar_dscale;
          factor_uscale = polar_uscale;
        } else {
          factor_pscale = special_polar_pscale[sbmask15(jextra)];
          factor_dscale = factor_uscale = 1.0;
        }

      } else {
        factor_wscale = special_polar_wscale[sbmask15(jextra)];
        if (igroup == jgroup) {
          factor_dscale = factor_pscale = special_polar_piscale[sbmask15(jextra)];
          factor_uscale = polar_uscale;
        } else {
          factor_dscale = factor_pscale = special_polar_pscale[sbmask15(jextra)];
          factor_uscale = 1.0;
        }
      }

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
      ukx = uind[j][0];
      uky = uind[j][1];
      ukz = uind[j][2];
      ukxp = uinp[j][0];
      ukyp = uinp[j][1];
      ukzp = uinp[j][2];

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
      uir = uix*xr + uiy*yr + uiz*zr;
      uirp = uixp*xr + uiyp*yr + uizp*zr;
      ukr = ukx*xr + uky*yr + ukz*zr;
      ukrp = ukxp*xr + ukyp*yr + ukzp*zr;

      // get reciprocal distance terms for this interaction

      rr1 = felec / r;
      rr3 = rr1 / r2;
      rr5 = 3.0 * rr3 / r2;
      rr7 = 5.0 * rr5 / r2;
      rr9 = 7.0 * rr7 / r2;

      // calculate the real space Ewald error function terms

      ralpha = aewald * r;
      bn[0] = erfc(ralpha) / r;
      alsq2 = 2.0 * aewald*aewald;
      alsq2n = 0.0;
      if (aewald > 0.0) alsq2n = 1.0 / (MY_PIS*aewald);
      exp2a = exp(-ralpha*ralpha);

      for (m = 1; m <= 4; m++) {
        bfac = (double) (m+m-1);
        alsq2n = alsq2 * alsq2n;
        bn[m] = (bfac*bn[m-1]+alsq2n*exp2a) / r2;
      }
      for (m = 0; m < 5; m++) bn[m] *= felec;

      // apply Thole polarization damping to scale factors

      sc3 = 1.0;
      sc5 = 1.0;
      sc7 = 1.0;
      for (k = 0; k < 3; k++) {
        rc3[k] = 0.0;
        rc5[k] = 0.0;
        rc7[k] = 0.0;
      }

      // apply Thole polarization damping to scale factors

      if (amoeba) {
        damp = pdi * pdamp[jtype];
        if (damp != 0.0) {
          pgamma = MIN(pti,thole[jtype]);
          damp = pgamma * cube(r/damp);
          if (damp < 50.0) {
            expdamp = exp(-damp);
            sc3 = 1.0 - expdamp;
            sc5 = 1.0 - (1.0+damp)*expdamp;
            sc7 = 1.0 - (1.0+damp+0.6*damp*damp) * expdamp;
            temp3 = 3.0 * damp * expdamp / r2;
            temp5 = damp;
            temp7 = -0.2 + 0.6*damp;
            rc3[0] = xr * temp3;
            rc3[1] = yr * temp3;
            rc3[2] = zr * temp3;
            rc5[0] = rc3[0] * temp5;
            rc5[1] = rc3[1] * temp5;
            rc5[2] = rc3[2] * temp5;
            rc7[0] = rc5[0] * temp7;
            rc7[1] = rc5[1] * temp7;
            rc7[2] = rc5[2] * temp7;
          }

          psc3 = 1.0 - sc3*factor_pscale;
          psc5 = 1.0 - sc5*factor_pscale;
          psc7 = 1.0 - sc7*factor_pscale;
          dsc3 = 1.0 - sc3*factor_dscale;
          dsc5 = 1.0 - sc5*factor_dscale;
          dsc7 = 1.0 - sc7*factor_dscale;
          usc3 = 1.0 - sc3*factor_uscale;
          usc5 = 1.0 - sc5*factor_uscale;
          psr3 = bn[1] - psc3*rr3;
          psr5 = bn[2] - psc5*rr5;
          psr7 = bn[3] - psc7*rr7;
          dsr3 = bn[1] - dsc3*rr3;
          dsr5 = bn[2] - dsc5*rr5;
          dsr7 = bn[3] - dsc7*rr7;
          usr5 = bn[2] - usc5*rr5;
          for (k = 0; k < 3; k++) {
            prc3[k] = rc3[k] * factor_pscale;
            prc5[k] = rc5[k] * factor_pscale;
            prc7[k] = rc7[k] * factor_pscale;
            drc3[k] = rc3[k] * factor_dscale;
            drc5[k] = rc5[k] * factor_dscale;
            drc7[k] = rc7[k] * factor_dscale;
            urc3[k] = rc3[k] * factor_uscale;
            urc5[k] = rc5[k] * factor_uscale;
          }
        } else {
          // avoid uninitialized data access when damp == 0.0
          psc3 = psc5 = psc7 = dsc3 = dsc5 = dsc7 = usc3 = usc5 = 0.0;
          psr3 = psr5 = psr7 = dsr3 = dsr5 = dsr7 = usr5 = 0.0;
          prc3[0] = prc3[1] = prc3[2] = 0.0;
          drc3[0] = drc3[1] = drc3[2] = 0.0;
          prc5[0] = prc5[1] = prc5[2] = 0.0;
          drc5[0] = drc5[1] = drc5[2] = 0.0;
          prc7[0] = prc7[1] = prc7[2] = 0.0;
          drc7[0] = drc7[1] = drc7[2] = 0.0;
          urc3[0] = urc3[1] = urc3[2] = 0.0;
          urc5[0] = urc5[1] = urc5[2] = 0.0;
        }

      // apply charge penetration damping to scale factors

      } else {
        corek = pcore[jclass];
        alphak = palpha[jclass];
        valk = pval[j];
        damppole(r,9,alphai,alphak,dmpi,dmpk,dmpik);
        rr3core = bn[1] - (1.0-factor_dscale)*rr3;
        rr5core = bn[2] - (1.0-factor_dscale)*rr5;
        rr3i = bn[1] - (1.0-factor_dscale*dmpi[2])*rr3;
        rr5i = bn[2] - (1.0-factor_dscale*dmpi[4])*rr5;
        rr7i = bn[3] - (1.0-factor_dscale*dmpi[6])*rr7;
        rr9i = bn[4] - (1.0-factor_dscale*dmpi[8])*rr9;
        rr3k = bn[1] - (1.0-factor_dscale*dmpk[2])*rr3;
        rr5k = bn[2] - (1.0-factor_dscale*dmpk[4])*rr5;
        rr7k = bn[3] - (1.0-factor_dscale*dmpk[6])*rr7;
        rr9k = bn[4] - (1.0-factor_dscale*dmpk[8])*rr9;
        rr5ik = bn[2] - (1.0-factor_wscale*dmpik[4])*rr5;
        rr7ik = bn[3] - (1.0-factor_wscale*dmpik[6])*rr7;
      }

      // get the induced dipole field used for dipole torques

      if (amoeba) {
        tix3 = psr3*ukx + dsr3*ukxp;
        tiy3 = psr3*uky + dsr3*ukyp;
        tiz3 = psr3*ukz + dsr3*ukzp;
        tkx3 = psr3*uix + dsr3*uixp;
        tky3 = psr3*uiy + dsr3*uiyp;
        tkz3 = psr3*uiz + dsr3*uizp;
        tuir = -psr5*ukr - dsr5*ukrp;
        tukr = -psr5*uir - dsr5*uirp;
      } else {
        tix3 = 2.0*rr3i*ukx;
        tiy3 = 2.0*rr3i*uky;
        tiz3 = 2.0*rr3i*ukz;
        tkx3 = 2.0*rr3k*uix;
        tky3 = 2.0*rr3k*uiy;
        tkz3 = 2.0*rr3k*uiz;
        tuir = -2.0*rr5i*ukr;
        tukr = -2.0*rr5k*uir;
      }

      ufld[i][0] += tix3 + xr*tuir;
      ufld[i][1] += tiy3 + yr*tuir;
      ufld[i][2] += tiz3 + zr*tuir;
      ufld[j][0] += tkx3 + xr*tukr;
      ufld[j][1] += tky3 + yr*tukr;
      ufld[j][2] += tkz3 + zr*tukr;

      // get induced dipole field gradient used for quadrupole torques

      if (amoeba) {
        tix5 = 2.0 * (psr5*ukx+dsr5*ukxp);
        tiy5 = 2.0 * (psr5*uky+dsr5*ukyp);
        tiz5 = 2.0 * (psr5*ukz+dsr5*ukzp);
        tkx5 = 2.0 * (psr5*uix+dsr5*uixp);
        tky5 = 2.0 * (psr5*uiy+dsr5*uiyp);
        tkz5 = 2.0 * (psr5*uiz+dsr5*uizp);
        tuir = -psr7*ukr - dsr7*ukrp;
        tukr = -psr7*uir - dsr7*uirp;

      } else {
        tix5 = 4.0 * (rr5i*ukx);
        tiy5 = 4.0 * (rr5i*uky);
        tiz5 = 4.0 * (rr5i*ukz);
        tkx5 = 4.0 * (rr5k*uix);
        tky5 = 4.0 * (rr5k*uiy);
        tkz5 = 4.0 * (rr5k*uiz);
        tuir = -2.0*rr7i*ukr;
        tukr = -2.0*rr7k*uir;
      }

      dufld[i][0] += xr*tix5 + xr*xr*tuir;
      dufld[i][1] += xr*tiy5 + yr*tix5 + 2.0*xr*yr*tuir;
      dufld[i][2] += yr*tiy5 + yr*yr*tuir;
      dufld[i][3] += xr*tiz5 + zr*tix5 + 2.0*xr*zr*tuir;
      dufld[i][4] += yr*tiz5 + zr*tiy5 + 2.0*yr*zr*tuir;
      dufld[i][5] += zr*tiz5 + zr*zr*tuir;

      dufld[j][0] -= xr*tkx5 + xr*xr*tukr;
      dufld[j][1] -= xr*tky5 + yr*tkx5 + 2.0*xr*yr*tukr;
      dufld[j][2] -= yr*tky5 + yr*yr*tukr;
      dufld[j][3] -= xr*tkz5 + zr*tkx5 + 2.0*xr*zr*tukr;
      dufld[j][4] -= yr*tkz5 + zr*tky5 + 2.0*yr*zr*tukr;
      dufld[j][5] -= zr*tkz5 + zr*zr*tukr;

      // get the dEd/dR terms used for direct polarization force

      if (amoeba) {
        term1 = bn[2] - dsc3*rr5;
        term2 = bn[3] - dsc5*rr7;
        term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3[0];
        term4 = rr3*drc3[0] - term1*xr - dsr5*xr;
        term5 = term2*xr*xr - dsr5 - rr5*xr*drc5[0];
        term6 = (bn[4]-dsc7*rr9)*xr*xr - bn[3] - rr7*xr*drc7[0];
        term7 = rr5*drc5[0] - 2.0*bn[3]*xr + (dsc5+1.5*dsc7)*rr7*xr;
        tixx = ci*term3 + dix*term4 + dir*term5 +
          2.0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7 + 2.0*qix*term7 + qir*term6;
        tkxx = ck*term3 - dkx*term4 - dkr*term5 +
          2.0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7 + 2.0*qkx*term7 + qkr*term6;
        term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3[1];
        term4 = rr3*drc3[1] - term1*yr - dsr5*yr;
        term5 = term2*yr*yr - dsr5 - rr5*yr*drc5[1];
        term6 = (bn[4]-dsc7*rr9)*yr*yr - bn[3] - rr7*yr*drc7[1];
        term7 = rr5*drc5[1] - 2.0*bn[3]*yr + (dsc5+1.5*dsc7)*rr7*yr;
        tiyy = ci*term3 + diy*term4 + dir*term5 +
          2.0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7 + 2.0*qiy*term7 + qir*term6;
        tkyy = ck*term3 - dky*term4 - dkr*term5 +
          2.0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7 + 2.0*qky*term7 + qkr*term6;
        term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3[2];
        term4 = rr3*drc3[2] - term1*zr - dsr5*zr;
        term5 = term2*zr*zr - dsr5 - rr5*zr*drc5[2];
        term6 = (bn[4]-dsc7*rr9)*zr*zr - bn[3] - rr7*zr*drc7[2];
        term7 = rr5*drc5[2] - 2.0*bn[3]*zr + (dsc5+1.5*dsc7)*rr7*zr;
        tizz = ci*term3 + diz*term4 + dir*term5 +
          2.0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7 + 2.0*qiz*term7 + qir*term6;
        tkzz = ck*term3 - dkz*term4 - dkr*term5 +
          2.0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7 + 2.0*qkz*term7 + qkr*term6;
        term3 = term1*xr*yr - rr3*yr*drc3[0];
        term4 = rr3*drc3[0] - term1*xr;
        term5 = term2*xr*yr - rr5*yr*drc5[0];
        term6 = (bn[4]-dsc7*rr9)*xr*yr - rr7*yr*drc7[0];
        term7 = rr5*drc5[0] - term2*xr;
        tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5 +
          2.0*dsr5*qixy - 2.0*dsr7*yr*qix + 2.0*qiy*term7 + qir*term6;
        tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5 +
          2.0*dsr5*qkxy - 2.0*dsr7*yr*qkx + 2.0*qky*term7 + qkr*term6;
        term3 = term1*xr*zr - rr3*zr*drc3[0];
        term5 = term2*xr*zr - rr5*zr*drc5[0];
        term6 = (bn[4]-dsc7*rr9)*xr*zr - rr7*zr*drc7[0];
        tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5 +
          2.0*dsr5*qixz - 2.0*dsr7*zr*qix + 2.0*qiz*term7 + qir*term6;
        tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5 +
          2.0*dsr5*qkxz - 2.0*dsr7*zr*qkx + 2.0*qkz*term7 + qkr*term6;
        term3 = term1*yr*zr - rr3*zr*drc3[1];
        term4 = rr3*drc3[1] - term1*yr;
        term5 = term2*yr*zr - rr5*zr*drc5[1];
        term6 = (bn[4]-dsc7*rr9)*yr*zr - rr7*zr*drc7[1];
        term7 = rr5*drc5[1] - term2*yr;
        tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5 +
          2.0*dsr5*qiyz - 2.0*dsr7*zr*qiy + 2.0*qiz*term7 + qir*term6;
        tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5 +
          2.0*dsr5*qkyz - 2.0*dsr7*zr*qky + 2.0*qkz*term7 + qkr*term6;
        depx = tixx*ukxp + tixy*ukyp + tixz*ukzp - tkxx*uixp - tkxy*uiyp - tkxz*uizp;
        depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp - tkxy*uixp - tkyy*uiyp - tkyz*uizp;
        depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp - tkxz*uixp - tkyz*uiyp - tkzz*uizp;
        frcx = depx;
        frcy = depy;
        frcz = depz;

        // get the dEp/dR terms used for direct polarization force

        term1 = bn[2] - psc3*rr5;
        term2 = bn[3] - psc5*rr7;
        term3 = -psr3 + term1*xr*xr - rr3*xr*prc3[0];
        term4 = rr3*prc3[0] - term1*xr - psr5*xr;
        term5 = term2*xr*xr - psr5 - rr5*xr*prc5[0];
        term6 = (bn[4]-psc7*rr9)*xr*xr - bn[3] - rr7*xr*prc7[0];
        term7 = rr5*prc5[0] - 2.0*bn[3]*xr + (psc5+1.5*psc7)*rr7*xr;
        tixx = ci*term3 + dix*term4 + dir*term5 +
          2.0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7 + 2.0*qix*term7 + qir*term6;
        tkxx = ck*term3 - dkx*term4 - dkr*term5 +
          2.0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7 + 2.0*qkx*term7 + qkr*term6;
        term3 = -psr3 + term1*yr*yr - rr3*yr*prc3[1];
        term4 = rr3*prc3[1] - term1*yr - psr5*yr;
        term5 = term2*yr*yr - psr5 - rr5*yr*prc5[1];
        term6 = (bn[4]-psc7*rr9)*yr*yr - bn[3] - rr7*yr*prc7[1];
        term7 = rr5*prc5[1] - 2.0*bn[3]*yr + (psc5+1.5*psc7)*rr7*yr;
        tiyy = ci*term3 + diy*term4 + dir*term5 +
          2.0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7 + 2.0*qiy*term7 + qir*term6;
        tkyy = ck*term3 - dky*term4 - dkr*term5 +
          2.0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7 + 2.0*qky*term7 + qkr*term6;
        term3 = -psr3 + term1*zr*zr - rr3*zr*prc3[2];
        term4 = rr3*prc3[2] - term1*zr - psr5*zr;
        term5 = term2*zr*zr - psr5 - rr5*zr*prc5[2];
        term6 = (bn[4]-psc7*rr9)*zr*zr - bn[3] - rr7*zr*prc7[2];
        term7 = rr5*prc5[2] - 2.0*bn[3]*zr + (psc5+1.5*psc7)*rr7*zr;
        tizz = ci*term3 + diz*term4 + dir*term5 +
          2.0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7 + 2.0*qiz*term7 + qir*term6;
        tkzz = ck*term3 - dkz*term4 - dkr*term5 +
          2.0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7 + 2.0*qkz*term7 + qkr*term6;
        term3 = term1*xr*yr - rr3*yr*prc3[0];
        term4 = rr3*prc3[0] - term1*xr;
        term5 = term2*xr*yr - rr5*yr*prc5[0];
        term6 = (bn[4]-psc7*rr9)*xr*yr - rr7*yr*prc7[0];
        term7 = rr5*prc5[0] - term2*xr;
        tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5 +
          2.0*psr5*qixy - 2.0*psr7*yr*qix + 2.0*qiy*term7 + qir*term6;
        tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5 +
          2.0*psr5*qkxy - 2.0*psr7*yr*qkx + 2.0*qky*term7 + qkr*term6;
        term3 = term1*xr*zr - rr3*zr*prc3[0];
        term5 = term2*xr*zr - rr5*zr*prc5[0];
        term6 = (bn[4]-psc7*rr9)*xr*zr - rr7*zr*prc7[0];
        tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5 +
          2.0*psr5*qixz - 2.0*psr7*zr*qix + 2.0*qiz*term7 + qir*term6;
        tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5 +
          2.0*psr5*qkxz - 2.0*psr7*zr*qkx + 2.0*qkz*term7 + qkr*term6;
        term3 = term1*yr*zr - rr3*zr*prc3[1];
        term4 = rr3*prc3[1] - term1*yr;
        term5 = term2*yr*zr - rr5*zr*prc5[1];
        term6 = (bn[4]-psc7*rr9)*yr*zr - rr7*zr*prc7[1];
        term7 = rr5*prc5[1] - term2*yr;
        tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5 +
          2.0*psr5*qiyz - 2.0*psr7*zr*qiy + 2.0*qiz*term7 + qir*term6;
        tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5 +
          2.0*psr5*qkyz - 2.0*psr7*zr*qky + 2.0*qkz*term7 + qkr*term6;
        depx = tixx*ukx + tixy*uky + tixz*ukz - tkxx*uix - tkxy*uiy - tkxz*uiz;
        depy = tixy*ukx + tiyy*uky + tiyz*ukz - tkxy*uix - tkyy*uiy - tkyz*uiz;
        depz = tixz*ukx + tiyz*uky + tizz*ukz - tkxz*uix - tkyz*uiy - tkzz*uiz;
        frcx = frcx + depx;
        frcy = frcy + depy;
        frcz = frcz + depz;

      // get the field gradient for direct polarization force

      } else {
        term1i = rr3i - rr5i*xr*xr;
        term1core = rr3core - rr5core*xr*xr;
        term2i = 2.0*rr5i*xr ;
        term3i = rr7i*xr*xr - rr5i;
        term4i = 2.0*rr5i;
        term5i = 5.0*rr7i*xr;
        term6i = rr9i*xr*xr;
        term1k = rr3k - rr5k*xr*xr;
        term2k = 2.0*rr5k*xr;
        term3k = rr7k*xr*xr - rr5k;
        term4k = 2.0*rr5k;
        term5k = 5.0*rr7k*xr;
        term6k = rr9k*xr*xr;
        tixx = vali*term1i + corei*term1core + dix*term2i - dir*term3i -
          qixx*term4i + qix*term5i - qir*term6i + (qiy*yr+qiz*zr)*rr7i;
        tkxx = valk*term1k + corek*term1core - dkx*term2k + dkr*term3k -
          qkxx*term4k + qkx*term5k - qkr*term6k + (qky*yr+qkz*zr)*rr7k;
        term1i = rr3i - rr5i*yr*yr;
        term1core = rr3core - rr5core*yr*yr;
        term2i = 2.0*rr5i*yr;
        term3i = rr7i*yr*yr - rr5i;
        term4i = 2.0*rr5i;
        term5i = 5.0*rr7i*yr;
        term6i = rr9i*yr*yr;
        term1k = rr3k - rr5k*yr*yr;
        term2k = 2.0*rr5k*yr;
        term3k = rr7k*yr*yr - rr5k;
        term4k = 2.0*rr5k;
        term5k = 5.0*rr7k*yr;
        term6k = rr9k*yr*yr;
        tiyy = vali*term1i + corei*term1core + diy*term2i - dir*term3i -
          qiyy*term4i + qiy*term5i - qir*term6i + (qix*xr+qiz*zr)*rr7i;
        tkyy = valk*term1k + corek*term1core - dky*term2k + dkr*term3k -
          qkyy*term4k + qky*term5k - qkr*term6k + (qkx*xr+qkz*zr)*rr7k;
        term1i = rr3i - rr5i*zr*zr;
        term1core = rr3core - rr5core*zr*zr;
        term2i = 2.0*rr5i*zr;
        term3i = rr7i*zr*zr - rr5i;
        term4i = 2.0*rr5i;
        term5i = 5.0*rr7i*zr;
        term6i = rr9i*zr*zr;
        term1k = rr3k - rr5k*zr*zr;
        term2k = 2.0*rr5k*zr;
        term3k = rr7k*zr*zr - rr5k;
        term4k = 2.0*rr5k;
        term5k = 5.0*rr7k*zr;
        term6k = rr9k*zr*zr;
        tizz = vali*term1i + corei*term1core + diz*term2i - dir*term3i -
          qizz*term4i + qiz*term5i - qir*term6i + (qix*xr+qiy*yr)*rr7i;
        tkzz = valk*term1k + corek*term1core - dkz*term2k + dkr*term3k -
          qkzz*term4k + qkz*term5k - qkr*term6k + (qkx*xr+qky*yr)*rr7k;
        term2i = rr5i*xr ;
        term1i = yr * term2i;
        term1core = rr5core*xr*yr;
        term3i = rr5i*yr;
        term4i = yr * (rr7i*xr);
        term5i = 2.0*rr5i;
        term6i = 2.0*rr7i*xr;
        term7i = 2.0*rr7i*yr;
        term8i = yr*rr9i*xr;
        term2k = rr5k*xr;
        term1k = yr * term2k;
        term3k = rr5k*yr;
        term4k = yr * (rr7k*xr);
        term5k = 2.0*rr5k;
        term6k = 2.0*rr7k*xr;
        term7k = 2.0*rr7k*yr;
        term8k = yr*rr9k*xr;
        tixy = -vali*term1i - corei*term1core + diy*term2i + dix*term3i -
          dir*term4i - qixy*term5i + qiy*term6i + qix*term7i - qir*term8i;
        tkxy = -valk*term1k - corek*term1core - dky*term2k - dkx*term3k +
          dkr*term4k - qkxy*term5k + qky*term6k + qkx*term7k - qkr*term8k;
        term2i = rr5i*xr;
        term1i = zr * term2i;
        term1core = rr5core*xr*zr;
        term3i = rr5i*zr;
        term4i = zr * (rr7i*xr);
        term5i = 2.0*rr5i;
        term6i = 2.0*rr7i*xr;
        term7i = 2.0*rr7i*zr;
        term8i = zr*rr9i*xr;
        term2k = rr5k*xr;
        term1k = zr * term2k;
        term3k = rr5k*zr;
        term4k = zr * (rr7k*xr);
        term5k = 2.0*rr5k;
        term6k = 2.0*rr7k*xr;
        term7k = 2.0*rr7k*zr;
        term8k = zr*rr9k*xr;
        tixz = -vali*term1i - corei*term1core + diz*term2i + dix*term3i -
          dir*term4i - qixz*term5i + qiz*term6i + qix*term7i - qir*term8i;
        tkxz = -valk*term1k - corek*term1core - dkz*term2k - dkx*term3k +
          dkr*term4k - qkxz*term5k + qkz*term6k + qkx*term7k - qkr*term8k;
        term2i = rr5i*yr;
        term1i = zr * term2i;
        term1core = rr5core*yr*zr;
        term3i = rr5i*zr;
        term4i = zr * (rr7i*yr);
        term5i = 2.0*rr5i;
        term6i = 2.0*rr7i*yr;
        term7i = 2.0*rr7i*zr;
        term8i = zr*rr9i*yr;
        term2k = rr5k*yr;
        term1k = zr * term2k;
        term3k = rr5k*zr;
        term4k = zr * (rr7k*yr);
        term5k = 2.0*rr5k;
        term6k = 2.0*rr7k*yr;
        term7k = 2.0*rr7k*zr;
        term8k = zr*rr9k*yr;
        tiyz = -vali*term1i - corei*term1core + diz*term2i + diy*term3i -
          dir*term4i - qiyz*term5i + qiz*term6i + qiy*term7i - qir*term8i;
        tkyz = -valk*term1k - corek*term1core - dkz*term2k - dky*term3k +
          dkr*term4k - qkyz*term5k + qkz*term6k + qky*term7k - qkr*term8k;
        depx = tixx*ukx + tixy*uky + tixz*ukz - tkxx*uix - tkxy*uiy - tkxz*uiz;
        depy = tixy*ukx + tiyy*uky + tiyz*ukz - tkxy*uix - tkyy*uiy - tkyz*uiz;
        depz = tixz*ukx + tiyz*uky + tizz*ukz - tkxz*uix - tkyz*uiy - tkzz*uiz;
        frcx = -2.0 * depx;
        frcy = -2.0 * depy;
        frcz = -2.0 * depz;
      }

      // get the dtau/dr terms used for mutual polarization force

      if (poltyp == MUTUAL && amoeba) {
        term1 = bn[2] - usc3*rr5;
        term2 = bn[3] - usc5*rr7;
        term3 = usr5 + term1;
        term4 = rr3 * factor_uscale;
        term5 = -xr*term3 + rc3[0]*term4;
        term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5[0];
        tixx = uix*term5 + uir*term6;
        tkxx = ukx*term5 + ukr*term6;
        term5 = -yr*term3 + rc3[1]*term4;
        term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5[1];
        tiyy = uiy*term5 + uir*term6;
        tkyy = uky*term5 + ukr*term6;
        term5 = -zr*term3 + rc3[2]*term4;
        term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5[2];
        tizz = uiz*term5 + uir*term6;
        tkzz = ukz*term5 + ukr*term6;
        term4 = -usr5 * yr;
        term5 = -xr*term1 + rr3*urc3[0];
        term6 = xr*yr*term2 - rr5*yr*urc5[0];
        tixy = uix*term4 + uiy*term5 + uir*term6;
        tkxy = ukx*term4 + uky*term5 + ukr*term6;
        term4 = -usr5 * zr;
        term6 = xr*zr*term2 - rr5*zr*urc5[0];
        tixz = uix*term4 + uiz*term5 + uir*term6;
        tkxz = ukx*term4 + ukz*term5 + ukr*term6;
        term5 = -yr*term1 + rr3*urc3[1];
        term6 = yr*zr*term2 - rr5*zr*urc5[1];
        tiyz = uiy*term4 + uiz*term5 + uir*term6;
        tkyz = uky*term4 + ukz*term5 + ukr*term6;
        depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
          + tkxx*uixp + tkxy*uiyp + tkxz*uizp;
        depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
          + tkxy*uixp + tkyy*uiyp + tkyz*uizp;
        depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
          + tkxz*uixp + tkyz*uiyp + tkzz*uizp;
        frcx = frcx + depx;
        frcy = frcy + depy;
        frcz = frcz + depz;

      // get the dtau/dr terms used for mutual polarization force

      } else if (poltyp == MUTUAL && !amoeba) {
        term1 = 2.0 * rr5ik;
        term2 = term1*xr;
        term3 = rr5ik - rr7ik*xr*xr;
        tixx = uix*term2 + uir*term3;
        tkxx = ukx*term2 + ukr*term3;
        term2 = term1*yr;
        term3 = rr5ik - rr7ik*yr*yr;
        tiyy = uiy*term2 + uir*term3;
        tkyy = uky*term2 + ukr*term3;
        term2 = term1*zr;
        term3 = rr5ik - rr7ik*zr*zr;
        tizz = uiz*term2 + uir*term3;
        tkzz = ukz*term2 + ukr*term3;
        term1 = rr5ik*yr;
        term2 = rr5ik*xr;
        term3 = yr * (rr7ik*xr);
        tixy = uix*term1 + uiy*term2 - uir*term3;
        tkxy = ukx*term1 + uky*term2 - ukr*term3;
        term1 = rr5ik * zr;
        term3 = zr * (rr7ik*xr);
        tixz = uix*term1 + uiz*term2 - uir*term3;
        tkxz = ukx*term1 + ukz*term2 - ukr*term3;
        term2 = rr5ik*yr;
        term3 = zr * (rr7ik*yr);
        tiyz = uiy*term1 + uiz*term2 - uir*term3;
        tkyz = uky*term1 + ukz*term2 - ukr*term3;
        depx = tixx*ukxp + tixy*ukyp + tixz*ukzp + tkxx*uixp + tkxy*uiyp + tkxz*uizp;
        depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp + tkxy*uixp + tkyy*uiyp + tkyz*uizp;
        depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp + tkxz*uixp + tkyz*uiyp + tkzz*uizp;
        frcx = frcx - depx;
        frcy = frcy - depy;
        frcz = frcz - depz;

      // get the dtau/dr terms used for OPT polarization force

      } else if (poltyp == OPT && amoeba) {
        for (k = 0; k < optorder; k++) {
          uirm = uopt[i][k][0]*xr + uopt[i][k][1]*yr + uopt[i][k][2]*zr;
          for (m = 0; m < optorder-k; m++) {
            ukrm = uopt[j][m][0]*xr + uopt[j][m][1]*yr + uopt[j][m][2]*zr;
            term1 = bn[2] - usc3*rr5;
            term2 = bn[3] - usc5*rr7;
            term3 = usr5 + term1;
            term4 = rr3 * factor_uscale;
            term5 = -xr*term3 + rc3[0]*term4;
            term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5[0];
            tixx = uopt[i][k][0]*term5 + uirm*term6;
            tkxx = uopt[j][m][0]*term5 + ukrm*term6;
            term5 = -yr*term3 + rc3[1]*term4;
            term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5[1];
            tiyy = uopt[i][k][1]*term5 + uirm*term6;
            tkyy = uopt[j][m][1]*term5 + ukrm*term6;
            term5 = -zr*term3 + rc3[2]*term4;
            term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5[2];
            tizz = uopt[i][k][2]*term5 + uirm*term6;
            tkzz = uopt[j][m][2]*term5 + ukrm*term6;
            term4 = -usr5 * yr;
            term5 = -xr*term1 + rr3*urc3[0];
            term6 = xr*yr*term2 - rr5*yr*urc5[0];
            tixy = uopt[i][k][0]*term4 + uopt[i][k][1]*term5 + uirm*term6;
            tkxy = uopt[j][m][0]*term4 + uopt[j][m][1]*term5 + ukrm*term6;
            term4 = -usr5 * zr;
            term6 = xr*zr*term2 - rr5*zr*urc5[0];
            tixz = uopt[i][k][0]*term4 + uopt[i][k][2]*term5 + uirm*term6;
            tkxz = uopt[j][m][0]*term4 + uopt[j][m][2]*term5 + ukrm*term6;
            term5 = -yr*term1 + rr3*urc3[1];
            term6 = yr*zr*term2 - rr5*zr*urc5[1];
            tiyz = uopt[i][k][1]*term4 + uopt[i][k][2]*term5 + uirm*term6;
            tkyz = uopt[j][m][1]*term4 + uopt[j][m][2]*term5 + ukrm*term6;
            depx = tixx*uoptp[j][m][0] + tkxx*uoptp[i][k][0] + tixy*uoptp[j][m][1] +
              tkxy*uoptp[i][k][1] + tixz*uoptp[j][m][2] + tkxz*uoptp[i][k][2];
            depy = tixy*uoptp[j][m][0] + tkxy*uoptp[i][k][0] + tiyy*uoptp[j][m][1] +
              tkyy*uoptp[i][k][1] + tiyz*uoptp[j][m][2] + tkyz*uoptp[i][k][2];
            depz = tixz*uoptp[j][m][0] + tkxz*uoptp[i][k][0] + tiyz*uoptp[j][m][1] +
              tkyz*uoptp[i][k][1] + tizz*uoptp[j][m][2] + tkzz*uoptp[i][k][2];
            frcx += copm[k+m+1]*depx;
            frcy += copm[k+m+1]*depy;
            frcz += copm[k+m+1]*depz;
          }
        }

      // get the dtau/dr terms used for OPT polarization force

      } else if (poltyp == OPT && !amoeba) {
        for (k = 0; k < optorder; k++) {
          uirm = uopt[i][k][0]*xr + uopt[i][k][1]*yr + uopt[i][k][2]*zr;
          for (m = 0; m < optorder-k; m++) {
            ukrm = uopt[j][m][0]*xr + uopt[j][m][1]*yr + uopt[j][m][2]*zr;
            term1 = 2.0 * rr5ik;
            term2 = term1*xr;
            term3 = rr5ik - rr7ik*xr*xr;
            tixx = uopt[i][k][0]*term2 + uirm*term3;
            tkxx = uopt[j][m][0]*term2 + ukrm*term3;
            term2 = term1*yr;
            term3 = rr5ik - rr7ik*yr*yr;
            tiyy = uopt[i][k][1]*term2 + uirm*term3;
            tkyy = uopt[j][m][1]*term2 + ukrm*term3;
            term2 = term1*zr;
            term3 = rr5ik - rr7ik*zr*zr;
            tizz = uopt[i][k][2]*term2 + uirm*term3;
            tkzz = uopt[j][m][2]*term2 + ukrm*term3;
            term1 = rr5ik*yr;
            term2 = rr5ik*xr;
            term3 = yr * (rr7ik*xr);
            tixy = uopt[i][k][0]*term1 + uopt[i][k][1]*term2 - uirm*term3;
            tkxy = uopt[j][m][0]*term1 + uopt[j][m][1]*term2 - ukrm*term3;
            term1 = rr5ik * zr;
            term3 = zr * (rr7ik*xr);
            tixz = uopt[i][k][0]*term1 + uopt[i][k][2]*term2 - uirm*term3;
            tkxz = uopt[j][m][0]*term1 + uopt[j][m][2]*term2 - ukrm*term3;
            term2 = rr5ik*yr;
            term3 = zr * (rr7ik*yr);
            tiyz = uopt[i][k][1]*term1 + uopt[i][k][2]*term2 - uirm*term3;
            tkyz = uopt[j][m][1]*term1 + uopt[j][m][2]*term2 - ukrm*term3;
            depx = tixx*uoptp[j][m][0] + tkxx*uoptp[i][k][0] + tixy*uoptp[j][m][1] +
              tkxy*uoptp[i][k][1] + tixz*uoptp[j][m][2] + tkxz*uoptp[i][k][2];
            depy = tixy*uoptp[j][m][0] + tkxy*uoptp[i][k][0] + tiyy*uoptp[j][m][1] +
              tkyy*uoptp[i][k][1] + tiyz*uoptp[j][m][2] + tkyz*uoptp[i][k][2];
            depz = tixz*uoptp[j][m][0] + tkxz*uoptp[i][k][0] + tiyz*uoptp[j][m][1] +
              tkyz*uoptp[i][k][1] + tizz*uoptp[j][m][2] + tkzz*uoptp[i][k][2];
            frcx -= copm[k+m+1]*depx;
            frcy -= copm[k+m+1]*depy;
            frcz -= copm[k+m+1]*depz;
          }
        }

      // get the dtau/dr terms used for TCG polarization force

      } else if (poltyp == TCG) {
#if 0
        // poltyp TCG not yet supported for AMOEBA/HIPPO
        for (m = 0; m < tcgnab; m++) {
          ukx = ubd[m][j][0];
          uky = ubd[m][j][1];
          ukz = ubd[m][j][2];
          ukxp = ubp[m][j][0];
          ukyp = ubp[m][j][1];
          ukzp = ubp[m][j][2];
          uirt = uax[m]*xr + uay[m]*yr + uaz[m]*zr;
          ukrt = ukx*xr + uky*yr + ukz*zr;
          term1 = bn[2] - usc3*rr5;
          term2 = bn[3] - usc5*rr7;
          term3 = usr5 + term1;
          term4 = rr3 * factor_uscale;
          term5 = -xr*term3 + rc3[0]*term4;
          term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5[0];
          tixx = uax[m]*term5 + uirt*term6;
          tkxx = ukx*term5 + ukrt*term6;
          term5 = -yr*term3 + rc3[1]*term4;
          term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5[1];
          tiyy = uay[m]*term5 + uirt*term6;
          tkyy = uky*term5 + ukrt*term6;
          term5 = -zr*term3 + rc3[2]*term4;
          term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5[2];
          tizz = uaz[m]*term5 + uirt*term6;
          tkzz = ukz*term5 + ukrt*term6;
          term4 = -usr5 * yr;
          term5 = -xr*term1 + rr3*urc3[0];
          term6 = xr*yr*term2 - rr5*yr*urc5[0];
          tixy = uax[m]*term4 + uay[m]*term5 + uirt*term6;
          tkxy = ukx*term4 + uky*term5 + ukrt*term6;
          term4 = -usr5 * zr;
          term6 = xr*zr*term2 - rr5*zr*urc5[0];
          tixz = uax[m]*term4 + uaz[m]*term5 + uirt*term6;
          tkxz = ukx*term4 + ukz*term5 + ukrt*term6;
          term5 = -yr*term1 + rr3*urc3[1];
          term6 = yr*zr*term2 - rr5*zr*urc5[1];
          tiyz = uay[m]*term4 + uaz[m]*term5 + uirt*term6;
          tkyz = uky*term4 + ukz*term5 + ukrt*term6;
          depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
            + tkxx*uaxp[m] + tkxy*uayp[m]
            + tkxz*uazp[m];
          depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
            + tkxy*uaxp[m] + tkyy*uayp[m]
            + tkyz*uazp[m];
          depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
            + tkxz*uaxp[m] + tkyz*uayp[m]
            + tkzz*uazp[m];
          frcx += depx;
          frcy += depy;
          frcz += depz;

          ukx = uad[m][j][0];
          uky = uad[m][j][1];
          ukz = uad[m][j][2];
          ukxp = uap[m][j][0];
          ukyp = uap[m][j][1];
          ukzp = uap[m][j][2];
          uirt = ubx[m]*xr + uby[m]*yr + ubz[m]*zr;
          ukrt = ukx*xr + uky*yr + ukz*zr;
          term1 = bn[2] - usc3*rr5;
          term2 = bn[3] - usc5*rr7;
          term3 = usr5 + term1;
          term4 = rr3 * factor_uscale;
          term5 = -xr*term3 + rc3[0]*term4;
          term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5[0];
          tixx = ubx[m]*term5 + uirt*term6;
          tkxx = ukx*term5 + ukrt*term6;
          term5 = -yr*term3 + rc3[1]*term4;
          term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5[1];
          tiyy = uby[m]*term5 + uirt*term6;
          tkyy = uky*term5 + ukrt*term6;
          term5 = -zr*term3 + rc3[2]*term4;
          term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5[2];
          tizz = ubz[m]*term5 + uirt*term6;
          tkzz = ukz*term5 + ukrt*term6;
          term4 = -usr5 * yr;
          term5 = -xr*term1 + rr3*urc3[0];
          term6 = xr*yr*term2 - rr5*yr*urc5[0];
          tixy = ubx[m]*term4 + uby[m]*term5 + uirt*term6;
          tkxy = ukx*term4 + uky*term5 + ukrt*term6;
          term4 = -usr5 * zr;
          term6 = xr*zr*term2 - rr5*zr*urc5[0];
          tixz = ubx[m]*term4 + ubz[m]*term5 + uirt*term6;
          tkxz = ukx*term4 + ukz*term5 + ukrt*term6;
          term5 = -yr*term1 + rr3*urc3[1];
          term6 = yr*zr*term2 - rr5*zr*urc5[1];
          tiyz = uby[m]*term4 + ubz[m]*term5 + uirt*term6;
          tkyz = uky*term4 + ukz*term5 + ukrt*term6;
          depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
            + tkxx*ubxp[m] + tkxy*ubyp[m]
            + tkxz*ubzp[m];
          depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
            + tkxy*ubxp[m] + tkyy*ubyp[m]
            + tkyz*ubzp[m];
          depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
            + tkxz*ubxp[m] + tkyz*ubyp[m]
            + tkzz*ubzp[m];
          frcx += depx;
          frcy += depy;
          frcz += depz;
#endif
      }

      // increment force-based gradient on the interaction sites

      f[i][0] += frcx;
      f[i][1] += frcy;
      f[i][2] += frcz;
      f[j][0] -= frcx;
      f[j][1] -= frcy;
      f[j][2] -= frcz;

      // increment the virial due to pairwise Cartesian forces

      if (vflag_global) {
        vxx = xr * frcx;
        vxy = 0.5 * (yr*frcx+xr*frcy);
        vxz = 0.5 * (zr*frcx+xr*frcz);
        vyy = yr * frcy;
        vyz = 0.5 * (zr*frcy+yr*frcz);
        vzz = zr * frcz;

        virpolar[0] -= vxx;
        virpolar[1] -= vyy;
        virpolar[2] -= vzz;
        virpolar[3] -= vxy;
        virpolar[4] -= vxz;
        virpolar[5] -= vyz;
      }
    }
  }

  // reverse comm to sum ufld,dufld from ghost atoms to owned atoms

  crstyle = UFLD;
  comm->reverse_comm(this);

  // torque is induced field and gradient cross permanent moments

  for (i = 0; i < nlocal; i++) {
    dix = rpole[i][1];
    diy = rpole[i][2];
    diz = rpole[i][3];
    qixx = rpole[i][4];
    qixy = rpole[i][5];
    qixz = rpole[i][6];
    qiyy = rpole[i][8];
    qiyz = rpole[i][9];
    qizz = rpole[i][12];
    tep[0] = diz*ufld[i][1] - diy*ufld[i][2] +
      qixz*dufld[i][1] - qixy*dufld[i][3] +
      2.0*qiyz*(dufld[i][2]-dufld[i][5]) + (qizz-qiyy)*dufld[i][4];
    tep[1] = dix*ufld[i][2] - diz*ufld[i][0] -
      qiyz*dufld[i][1] + qixy*dufld[i][4] +
      2.0*qixz*(dufld[i][5]-dufld[i][0]) + (qixx-qizz)*dufld[i][3];
    tep[2] = diy*ufld[i][0] - dix*ufld[i][1] +
      qiyz*dufld[i][3] - qixz*dufld[i][4] +
      2.0*qixy*(dufld[i][0]-dufld[i][2]) + (qiyy-qixx)*dufld[i][1];

    torque2force(i,tep,fix,fiy,fiz,f);

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

    virpolar[0] -= vxx;
    virpolar[1] -= vyy;
    virpolar[2] -= vzz;
    virpolar[3] -= vxy;
    virpolar[4] -= vxz;
    virpolar[5] -= vyz;
  }
}

/* ----------------------------------------------------------------------
   polar_kspace = KSpace portion of induced dipole polarization
   adapted from Tinker eprecip1() routine
 ------------------------------------------------------------------------- */

void PairAmoeba::polar_kspace()
{
  int i,j,k,m,n;
  int nhalf1,nhalf2,nhalf3;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  int j1,j2,j3;
  int ix,iy,iz;
  double eterm,felec;
  double r1,r2,r3;
  double h1,h2,h3;
  double f1,f2,f3;
  double xix,yix,zix;
  double xiy,yiy,ziy;
  double xiz,yiz,ziz;
  double vxx,vyy,vzz;
  double vxy,vxz,vyz;
  double volterm,denom;
  double hsq,expterm;
  double term,pterm;
  double vterm,struc2;
  double tep[3];
  double fix[3],fiy[3],fiz[3];
  double cphid[4],cphip[4];
  double a[3][3];    // indices not flipped vs Fortran

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
  pterm = square(MY_PI/aewald);
  volterm = MY_PI * volbox;

  // initialize variables required for the scalar summation

  felec = electric / am_dielectric;

  // remove scalar sum virial from prior multipole FFT
  // can only do this if multipoles were computed with same aeewald = apewald
  // else need to re-compute it via new long-range solve

  nfft1 = p_kspace->nx;
  nfft2 = p_kspace->ny;
  nfft3 = p_kspace->nz;
  bsorder = p_kspace->order;

  nhalf1 = (nfft1+1) / 2;
  nhalf2 = (nfft2+1) / 2;
  nhalf3 = (nfft3+1) / 2;

  nxlo = p_kspace->nxlo_fft;
  nxhi = p_kspace->nxhi_fft;
  nylo = p_kspace->nylo_fft;
  nyhi = p_kspace->nyhi_fft;
  nzlo = p_kspace->nzlo_fft;
  nzhi = p_kspace->nzhi_fft;

  // use previous results or compute new qfac and convolution

  if (aewald == aeewald) {
    vxx = -vmsave[0];
    vyy = -vmsave[1];
    vzz = -vmsave[2];
    vxy = -vmsave[3];
    vxz = -vmsave[4];
    vyz = -vmsave[5];

  } else {

    // setup stencil size and B-spline coefficients

    moduli();
    bspline_fill();

    // convert Cartesian multipoles to fractional coordinates

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values

    double ***gridpre = (double ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    double *gridfft = p_kspace->pre_convolution();

    // ---------------------
    // convolution operation
    // ---------------------

    // zero virial accumulation variables

    vxx = vyy = vzz = vxy = vxz = vyz = 0.0;

    // perform convolution on K-space points I own

    m = n = 0;
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
            if (hsq) expterm = exp(term) / denom;
            struc2 = gridfft[n]*gridfft[n] + gridfft[n+1]*gridfft[n+1];
            eterm = 0.5 * felec * expterm * struc2;
            vterm = (2.0/hsq) * (1.0-term) * eterm;
            vxx -= h1*h1*vterm - eterm;
            vyy -= h2*h2*vterm - eterm;
            vzz -= h3*h3*vterm - eterm;
            vxy -= h1*h2*vterm;
            vxz -= h1*h3*vterm;
            vyz -= h2*h3*vterm;
          }

          expterm = qfac[m++];
          gridfft[n] *= expterm;
          gridfft[n+1] *= expterm;
          n += 2;
        }
      }
    }

    // post-convolution operations including backward FFT
    // gridppost = my portion of 3d grid in brick decomp w/ ghost values

    double ***gridpost = (double ***) p_kspace->post_convolution();

    // get potential

    fphi_mpole(gridpost,fphi);

    for (i = 0; i < nlocal; i++) {
      for (k = 0; k < 20; k++)
        fphi[i][k] *= felec;
    }

    // convert field from fractional to Cartesian

    fphi_to_cphi(fphi,cphi);
  }

  // convert Cartesian induced dipoles to fractional coordinates

  for (i = 0; i < 3; i++) {
    a[0][i] = nfft1 * recip[0][i];
    a[1][i] = nfft2 * recip[1][i];
    a[2][i] = nfft3 * recip[2][i];
  }

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 3; j++) {
      fuind[i][j] = a[j][0]*uind[i][0] + a[j][1]*uind[i][1] + a[j][2]*uind[i][2];
      fuinp[i][j] = a[j][0]*uinp[i][0] + a[j][1]*uinp[i][1] + a[j][2]*uinp[i][2];
    }
  }

  // gridpre2 = my portion of 4d grid in brick decomp w/ ghost values

  double ****gridpre2 = (double ****) pc_kspace->zero();

  // map 2 values to grid

  grid_uind(fuind,fuinp,gridpre2);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomposition

  double *gridfft = pc_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  // use qfac values from above or from induce()

  m = n = 0;
  for (k = nzlo; k <= nzhi; k++) {
    for (j = nylo; j <= nyhi; j++) {
      for (i = nxlo; i <= nxhi; i++) {
        term = qfac[m++];
        gridfft[n] *= term;
        gridfft[n+1] *= term;
        n += 2;
      }
    }
  }

  // post-convolution operations including backward FFT
  // gridppost = my portion of 4d grid in brick decomp w/ ghost values

  double ****gridpost = (double ****) pc_kspace->post_convolution();

  // get potential

  fphi_uind(gridpost,fphid,fphip,fphidp);

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 10; j++) {
      fphid[i][j] = felec * fphid[i][j];
      fphip[i][j] = felec * fphip[i][j];
    }
    for (j = 0; j < 20; j++)
      fphidp[i][j] = felec * fphidp[i][j];
  }

  // increment the dipole polarization gradient contributions

  for (i = 0; i < nlocal; i++) {
    f1 = 0.0;
    f2 = 0.0;
    f3 = 0.0;
    for (k = 0; k < 3; k++) {
      j1 = deriv1[k+1];
      j2 = deriv2[k+1];
      j3 = deriv3[k+1];
      f1 += (fuind[i][k]+fuinp[i][k])*fphi[i][j1];
      f2 += (fuind[i][k]+fuinp[i][k])*fphi[i][j2];
      f3 += (fuind[i][k]+fuinp[i][k])*fphi[i][j3];
      if (poltyp == MUTUAL) {
        f1 += fuind[i][k]*fphip[i][j1] + fuinp[i][k]*fphid[i][j1];
        f2 += fuind[i][k]*fphip[i][j2] + fuinp[i][k]*fphid[i][j2];
        f3 += fuind[i][k]*fphip[i][j3] + fuinp[i][k]*fphid[i][j3];
      }
    }
    for (k = 0; k < 10; k++) {
      f1 += fmp[i][k]*fphidp[i][deriv1[k]];
      f2 += fmp[i][k]*fphidp[i][deriv2[k]];
      f3 += fmp[i][k]*fphidp[i][deriv3[k]];
    }
    f1 *= 0.5 * nfft1;
    f2 *= 0.5 * nfft2;
    f3 *= 0.5 * nfft3;
    h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;
    h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
    h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;
    f[i][0] -= h1;
    f[i][1] -= h2;
    f[i][2] -= h3;
  }

  // set the potential to be the induced dipole average

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < 10; j++)
      fphidp[i][j] *= 0.5;
  }

  fphi_to_cphi(fphidp,cphidp);

  // get the fractional to Cartesian transformation matrix

  frac_to_cart();

  // increment the dipole polarization virial contributions

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++) {
      cphid[j] = 0.0;
      cphip[j] = 0.0;
      for (k = 1; k < 4; k++) {
        cphid[j] += ftc[j][k]*fphid[i][k];
        cphip[j] += ftc[j][k]*fphip[i][k];
      }
    }

    vxx -= cmp[i][1]*cphidp[i][1] +
      0.5*((uind[i][0]+uinp[i][0])*cphi[i][1]);
    vyy -= cmp[i][2]*cphidp[i][2] +
      0.5*((uind[i][1]+uinp[i][1])*cphi[i][2]);
    vzz -= cmp[i][3]*cphidp[i][3] +
      0.5*((uind[i][2]+uinp[i][2])*cphi[i][3]);
    vxy -= 0.5*(cphidp[i][1]*cmp[i][2]+cphidp[i][2]*cmp[i][1]) +
      0.25*((uind[i][1]+uinp[i][1])*cphi[i][1] +
            (uind[i][0]+uinp[i][0])*cphi[i][2]);
    vyz -= 0.5*(cphidp[i][2]*cmp[i][3]+cphidp[i][3]*cmp[i][2]) +
      0.25*((uind[i][2]+uinp[i][2])*cphi[i][2] +
            (uind[i][1]+uinp[i][1])*cphi[i][3]);
    vxz -= 0.5*(cphidp[i][1]*cmp[i][3]+cphidp[i][3]*cmp[i][1]) +
      0.25*((uind[i][2]+uinp[i][2])*cphi[i][1] +
            (uind[i][0]+uinp[i][0])*cphi[i][3]);

    vxx -= 2.0*cmp[i][4]*cphidp[i][4] + cmp[i][7]*cphidp[i][7] +
      cmp[i][8]*cphidp[i][8];
    vyy -= 2.0*cmp[i][5]*cphidp[i][5] + cmp[i][7]*cphidp[i][7] +
      cmp[i][9]*cphidp[i][9];
    vzz -= 2.0*cmp[i][6]*cphidp[i][6] + cmp[i][8]*cphidp[i][8] +
      cmp[i][9]*cphidp[i][9];
    vxy -= (cmp[i][4]+cmp[i][5])*cphidp[i][7] +
      0.5*(cmp[i][7]*(cphidp[i][5]+cphidp[i][4]) +
           cmp[i][8]*cphidp[i][9]+cmp[i][9]*cphidp[i][8]);
    vyz -= (cmp[i][5]+cmp[i][6])*cphidp[i][9] +
      0.5*(cmp[i][9]*(cphidp[i][5]+cphidp[i][6]) +
           cmp[i][7]*cphidp[i][8]+cmp[i][8]*cphidp[i][7]);
    vxz -= (cmp[i][4]+cmp[i][6])*cphidp[i][8] +
      0.5*(cmp[i][8]*(cphidp[i][4]+cphidp[i][6]) +
           cmp[i][7]*cphidp[i][9]+cmp[i][9]*cphidp[i][7]);

    if (poltyp == MUTUAL) {
      vxx -= 0.5 * (cphid[1]*uinp[i][0]+cphip[1]*uind[i][0]);
      vyy -= 0.5 * (cphid[2]*uinp[i][1]+cphip[2]*uind[i][1]);
      vzz -= 0.5 * (cphid[3]*uinp[i][2]+cphip[3]*uind[i][2]);
      vxy -= 0.25 * (cphid[1]*uinp[i][1]+cphip[1]*uind[i][1] +
                     cphid[2]*uinp[i][0]+cphip[2]*uind[i][0]);
      vyz -= 0.25 * (cphid[2]*uinp[i][2]+cphip[2]*uind[i][2] +
                     cphid[3]*uinp[i][1]+cphip[3]*uind[i][1]);
      vxz -= 0.25 * (cphid[1]*uinp[i][2]+cphip[1]*uind[i][2] +
                     cphid[3]*uinp[i][0]+cphip[3]*uind[i][0]);
    }
  }


  // resolve site torques then increment forces and virial

  for (i = 0; i < nlocal; i++) {
    tep[0] = cmp[i][3]*cphidp[i][2] - cmp[i][2]*cphidp[i][3] +
      2.0*(cmp[i][6]-cmp[i][5])*cphidp[i][9] + cmp[i][8]*cphidp[i][7] +
      cmp[i][9]*cphidp[i][5]- cmp[i][7]*cphidp[i][8] - cmp[i][9]*cphidp[i][6];
    tep[1] = cmp[i][1]*cphidp[i][3] - cmp[i][3]*cphidp[i][1] +
      2.0*(cmp[i][4]-cmp[i][6])*cphidp[i][8] + cmp[i][7]*cphidp[i][9] +
      cmp[i][8]*cphidp[i][6] - cmp[i][8]*cphidp[i][4] - cmp[i][9]*cphidp[i][7];
    tep[2] = cmp[i][2]*cphidp[i][1] - cmp[i][1]*cphidp[i][2] +
      2.0*(cmp[i][5]-cmp[i][4])*cphidp[i][7] + cmp[i][7]*cphidp[i][4] +
      cmp[i][9]*cphidp[i][8] - cmp[i][7]*cphidp[i][5] - cmp[i][8]*cphidp[i][9];

    torque2force(i,tep,fix,fiy,fiz,f);

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
    vyy += yix*fix[1] + yiy*fiy[1] + yiz*fiz[1];
    vzz += zix*fix[2] + ziy*fiy[2] + ziz*fiz[2];
    vxy += 0.5*(yix*fix[0] + yiy*fiy[0] + yiz*fiz[0] +
                xix*fix[1] + xiy*fiy[1] + xiz*fiz[1]);
    vyz += 0.5*(zix*fix[1] + ziy*fiy[1] + ziz*fiz[1] +
                yix*fix[2] + yiy*fiy[2] + yiz*fiz[2]);
    vxz += 0.5*(zix*fix[0] + ziy*fiy[0] + ziz*fiz[0] +
                xix*fix[2] + xiy*fiy[2] + xiz*fiz[2]);
  }

  // account for dipole response terms in the OPT method

  if (poltyp == OPT) {
    for (i = 0; i < nlocal; i++) {
      for (k = 0; k < optorder; k++) {
        for (j = 1; j < 10; j++) {
          fphid[i][j] = felec * fopt[i][k][j];
          fphip[i][j] = felec * foptp[i][k][j];
        }

        for (m = 0; m < optorder-k; m++) {
          for (j = 0; j < 3; j++) {
            fuind[i][j] = a[0][j]*uopt[i][m][0] + a[1][j]*uopt[i][m][1] +
              a[2][j]*uopt[i][m][2];
            fuinp[i][j] = a[0][j]*uoptp[i][m][0] + a[1][j]*uoptp[i][m][1] +
              a[2][j]*uoptp[i][m][2];
          }

          f1 = 0.0;
          f2 = 0.0;
          f3 = 0.0;

          for (j = 0; j < 3; j++) {
            j1 = deriv1[j+1];
            j2 = deriv2[j+1];
            j3 = deriv3[j+1];
            f1 += fuind[i][j]*fphip[i][j1] + fuinp[i][j]*fphid[i][j1];
            f2 += fuind[i][j]*fphip[i][j2] + fuinp[i][j]*fphid[i][j2];
            f3 += fuind[i][j]*fphip[i][j3] + fuinp[i][j]*fphid[i][j3];
          }

          f1 *= 0.5 * nfft1;
          f2 *= 0.5 * nfft2;
          f3 *= 0.5 * nfft3;
          h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;
          h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
          h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;

          f[i][0] -= copm[k+m+1]*h1;
          f[i][1] -= copm[k+m+1]*h2;
          f[i][2] -= copm[k+m+1]*h3;

          for (j = 1; j < 4; j++) {
            cphid[j] = 0.0;
            cphip[j] = 0.0;
            for (j1 = 1; j1 < 4; j1++) {
              cphid[j] += ftc[j][j1]*fphid[i][j1];
              cphip[j] += ftc[j][j1]*fphip[i][j1];
            }
          }

          vxx -= 0.5*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][0] + cphip[1]*uopt[i][m][0]);
          vyy -= 0.5*copm[k+m+1] *
            (cphid[2]*uoptp[i][m][1]+ cphip[2]*uopt[i][m][1]);
          vzz -= 0.5*copm[k+m+1] *
            (cphid[3]*uoptp[i][m][2]+ cphip[3]*uopt[i][m][2]);
          vxy -= 0.25*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][1]+ cphip[1]*uopt[i][m][1]+
             cphid[2]*uoptp[i][m][0]+ cphip[2]*uopt[i][m][0]);
          vyz -= 0.25*copm[k+m+1] *
            (cphid[1]*uoptp[i][m][2]+ cphip[1]*uopt[i][m][2]+
             cphid[3]*uoptp[i][m][0]+ cphip[3]*uopt[i][m][0]);
          vxz -= 0.25*copm[k+m+1] *
            (cphid[2]*uoptp[i][m][2]+ cphip[2]*uopt[i][m][2]+
             cphid[3]*uoptp[i][m][1]+ cphip[3]*uopt[i][m][1]);
        }
      }
    }
  }

  // account for dipole response terms in the TCG method

  /*
  if (poltyp == TCG) {

    for (m = 0; m < tcgnab; m++) {
      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          fuind[i][j] = a[0][j]*uad[i][m][0] + a[1][j]*uad[i][m][1] +
            a[2][j]*uad[i][m][2];
          fuinp[i][j] = a[0][j]*ubp[i][m][0] + a[1][j]*ubp[i][m][1] +
            a[2][j]*ubp[i][m][2];
        }
      }

      grid_uind(fuind,fuinp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      for (k = 0; k < nfft3; k++) {
        for (j = 0; j < nfft2; j++) {
          for (i = 0; i < nfft1; i++) {
            term = qfac[k][j][i];
            qgrid[k][j][i][0] *= term;
            qgrid[k][j][i][1] *= term;
          }
        }
      }

      efft->compute(qgrid[0][0][0],qgrid[0][0][0],-1);
      fphi_uind(fphid,fphip,fphidp);

      for (i = 0; i < nlocal; i++) {
        for (j = 1; j < 10; j++) {
          fphid[i][j] *= felec;
          fphip[i][j] *= felec;
        }
      }

      for (i = 0; i < nlocal; i++) {
        f1 = 0.0;
        f2 = 0.0;
        f3 = 0.0;
        for (k = 0; k < 3; k++) {
          j1 = deriv1[k+1];
          j2 = deriv2[k+1];
          j3 = deriv3[k+1];
          f1 += fuind[i][k]*fphip[i][j1]+fuinp[i][k]*fphid[i][j1];
          f2 += fuind[i][k]*fphip[i][j2]+fuinp[i][k]*fphid[i][j2];
          f3 += fuind[i][k]*fphip[i][j3]+fuinp[i][k]*fphid[i][j3];
        }

        f1 *= 0.5 * nfft1;
        f2 *= 0.5 * nfft2;
        f3 *= 0.5 * nfft3;
        h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;
        h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
        h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;
        f[i][0] -= h1;
        f[i][1] -= h2;
        f[i][2] -= h3;

        for (j = 1; j < 4; j++) {
          cphid[j] = 0.0;
          cphip[j] = 0.0;
          for (k = 1; k < 4; k++) {
            cphid[j] += ftc[j][k]*fphid[i][k];
            cphip[j] += ftc[j][k]*fphip[i][k];
          }
        }

        vxx -= 0.5*(cphid[1]*ubp[i][m][0] + cphip[1]*uad[i][m][0]);
        vyy -= 0.5*(cphid[2]*ubp[i][m][1] + cphip[2]*uad[i][m][1]);
        vzz -= 0.5*(cphid[3]*ubp[i][m][2] + cphip[3]*uad[i][m][2]);

        vxy -= 0.25*(cphid[1]*ubp[i][m][1] + cphip[1]*uad[i][m][1] +
                        cphid[2]*ubp[i][m][0] + cphip[2]*uad[i][m][0]);
        vyz -= 0.25*(cphid[1]*ubp[i][m][2] + cphip[1]*uad[i][m][2] +
                        cphid[3]*ubp[i][m][0] + cphip[3]*uad[i][m][0]);
        vxz -= 0.25*(cphid[2]*ubp[i][m][2] + cphip[2]*uad[i][m][2] +
                        cphid[3]*ubp[i][m][1] + cphip[3]*uad[i][m][1]);
      }

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 3; j++) {
          fuind[i][j] = a[0][j]*ubd[i][m][0] + a[1][j]*ubd[i][m][1] +
            a[2][j]*ubd[i][m][2];
          fuinp[i][j] = a[0][j]*uap[i][m][0] + a[1][j]*uap[i][m][1] +
            a[2][j]*uap[i][m][2];
        }
      }

      grid_uind(fuind,fuinp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      for (k = 0; k < nfft3; k++) {
        for (j = 0; j < nfft2; j++) {
          for (i = 0; i < nfft1; i++) {
            term = qfac[k][j][i];
            qgrid[k][j][i][0] *= term;
            qgrid[k][j][i][1] *= term;
          }
        }
      }

      efft->compute(qgrid[0][0][0],qgrid[0][0][0],-1);
      fphi_uind(fphid,fphip,fphidp);

      for (i = 0; i < nlocal; i++) {
        for (j = 1; j < 10; j++) {
          fphid[i][j] *= felec;
          fphip[i][j] *= felec;
        }
      }

      for (i = 0; i < nlocal; i++) {
        f1 = 0.0;
        f2 = 0.0;
        f3 = 0.0;
        for (k = 0; k < 3; k++) {
          j1 = deriv1[k+1];
          j2 = deriv2[k+1];
          j3 = deriv3[k+1];
          f1 += fuind[i][k]*fphip[i][j1]+fuinp[i][k]*fphid[i][j1];
          f2 += fuind[i][k]*fphip[i][j2]+fuinp[i][k]*fphid[i][j2];
          f3 += fuind[i][k]*fphip[i][j3]+fuinp[i][k]*fphid[i][j3];
        }

        f1 *= 0.5 * nfft1;
        f2 *= 0.5 * nfft2;
        f3 *= 0.5 * nfft3;
        h1 = recip[0][0]*f1 + recip[0][1]*f2 + recip[0][2]*f3;  // matvec
        h2 = recip[1][0]*f1 + recip[1][1]*f2 + recip[1][2]*f3;
        h3 = recip[2][0]*f1 + recip[2][1]*f2 + recip[2][2]*f3;
        f[i][0] -= h1;
        f[i][1] -= h2;
        f[i][2] -= h3;

        for (j = 1; j < 4; j++) {
          cphid[j] = 0.0;
          cphip[j] = 0.0;
          for (k = 1; k < 4; k++) {
            cphid[j] += ftc[j][k]*fphid[i][k];
            cphip[j] += ftc[j][k]*fphip[i][k];
          }
        }

        vxx -= 0.5*(cphid[1]*uap[i][m][0] + cphip[1]*ubd[i][m][0]);
        vyy -= 0.5*(cphid[2]*uap[i][m][1] + cphip[2]*ubd[i][m][1]);
        vzz -= 0.5*(cphid[3]*uap[i][m][2] + cphip[3]*ubd[i][m][2]);
        vxy -= 0.25*(cphid[1]*uap[i][m][1] + cphip[1]*ubd[i][m][1] +
                     cphid[2]*uap[i][m][0] + cphip[2]*ubd[i][m][0]);
        vxz -= 0.25*(cphid[1]*uap[i][m][2] + cphip[1]*ubd[i][m][2] +
                     cphid[3]*uap[i][m][0] + cphip[3]*ubd[i][m][0]);
        vyz -= 0.25*(cphid[2]*uap[i][m][2] + cphip[2]*ubd[i][m][2] +
                     cphid[3]*uap[i][m][1] + cphip[3]*ubd[i][m][1]);
      }
    }
  }
  */

  // assign permanent and induced multipoles to the PME grid

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++)
      cmp[i][j] += uinp[i][j-1];
  }

  // convert Cartesian multipoles to fractional multipoles

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by zero()

  double ***gridpre = (double ***) p_kspace->zero();

  // map atoms to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

  gridfft = p_kspace->pre_convolution();

  // gridfft1 = copy of first FFT

  int nfft_owned = p_kspace->nfft_owned;
  memcpy(gridfft1,gridfft,2*nfft_owned*sizeof(FFT_SCALAR));

  // assign induced dipoles to the PME grid

  for (i = 0; i < nlocal; i++) {
    for (j = 1; j < 4; j++)
      cmp[i][j] += uind[i][j-1] - uinp[i][j-1];
  }

  // convert Cartesian multipoles to fractional multipoles

  cmp_to_fmp(cmp,fmp);

  // gridpre = my portion of 3d grid in brick decomp w/ ghost values
  // zeroed by zero()

  gridpre = (double ***) p_kspace->zero();

  // map atoms to grid

  grid_mpole(fmp,gridpre);

  // pre-convolution operations including forward FFT
  // gridfft1/2 = my portions of complex 3d grid in FFT decomp as 1d vectors

  double *gridfft2 = p_kspace->pre_convolution();

  // ---------------------
  // convolution operation
  // ---------------------

  m = n = 0;
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
          struc2 = gridfft1[n]*gridfft2[n] + gridfft1[n+1]*gridfft2[n+1];
          eterm = 0.5 * felec * expterm * struc2;
          vterm = (2.0/hsq) * (1.0-term) * eterm;
          vxx += h1*h1*vterm - eterm;
          vyy += h2*h2*vterm - eterm;
          vzz += h3*h3*vterm - eterm;
          vxy += h1*h2*vterm;
          vyz += h2*h3*vterm;
          vxz += h1*h3*vterm;
        }
        n += 2;
      }
    }
  }

  // assign only the induced dipoles to the PME grid
  // and perform the 3-D FFT forward transformation
  // NOTE: why is there no inverse FFT in this section?

  if (poltyp == DIRECT || poltyp == TCG) {

    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < 10; j++)
        cmp[i][j] = 0.0;
      for (j = 1; j < 4; j++)
        cmp[i][j] = uinp[i][j-1];
    }

    // convert Cartesian multipoles to fractional multipoles

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values
    // zeroed by zero()

    double ***gridpre = (double ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    double *gridfft = p_kspace->pre_convolution();

    // gridfft1 = copy of first FFT

    int nfft_owned = p_kspace->nfft_owned;
    memcpy(gridfft1,gridfft,2*nfft_owned*sizeof(double));

    // assign ??? to the PME grid

    for (i = 0; i < nlocal; i++) {
      for (j = 1; j < 4; j++)
        cmp[i][j] = uind[i][j-1];
    }

    // convert Cartesian multipoles to fractional multipoles

    cmp_to_fmp(cmp,fmp);

    // gridpre = my portion of 3d grid in brick decomp w/ ghost values

    gridpre = (double ***) p_kspace->zero();

    // map atoms to grid

    grid_mpole(fmp,gridpre);

    // pre-convolution operations including forward FFT
    // gridfft = my portion of complex 3d grid in FFT decomp as 1d vector

    double *gridfft2 = p_kspace->pre_convolution();

    // ---------------------
    // convolution operation
    // ---------------------

    m = n = 0;
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
            struc2 = gridfft1[n]*gridfft2[n] + gridfft1[n+1]*gridfft2[n+1];
            eterm = 0.5 * felec * expterm * struc2;
            vterm = (2.0/hsq) * (1.0-term) * eterm;
            vxx += h1*h1*vterm - eterm;
            vyy += h2*h2*vterm - eterm;
            vzz += h3*h3*vterm - eterm;
            vxy += h1*h2*vterm;
            vyz += h2*h3*vterm;
            vxz += h1*h3*vterm;
          }
          n += 2;
        }
      }
    }
  }

  // add back missing terms for the TCG polarization method;
  // first do the term for "UAD" dotted with "UBP"

  /*
  if (poltyp == TCG) {

    for (m = 0; m < tcgnab; m++) {
      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 10; j++)
          cmp[i][j] = 0.0;
        for (j = 1; j < 4; j++)
          cmp[i][j] = ubp[i][m][j-1];
      }

      cmp_to_fmp(cmp,fmp);
      grid_mpole(fmp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      for (k = 0; k < nfft3; k++) {
        for (j = 0; j < nfft2; j++) {
          for (i = 0; i < nfft1; i++) {
            qgrip[k][j][i][0] = qgrid[k][j][i][0];
            qgrip[k][j][i][1] = qgrid[k][j][i][1];
          }
        }
      }

      for (i = 0; i < nlocal; i++) {
        for (j = 1; j < 4; j++)
          cmp[i][j] = uad[i][m][j-1];
      }

      cmp_to_fmp(cmp,fmp);
      grid_mpole(fmp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      // make the scalar summation over reciprocal lattice
      // NOTE: this loop has to be distributed for parallel
      // NOTE: why does this one include m = 0 ?

      for (m = 1; m < ntot; m++) {
        k1 = m % nfft1;
        k2 = (m % nff) / nfft1;
        k3 = m/nff;
        r1 = (k1 >= nf1) ? k1-nfft1 : k1;
        r2 = (k2 >= nf2) ? k2-nfft2 : k2;
        r3 = (k3 >= nf3) ? k3-nfft3 : k3;
        h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;
        h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
        h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
        hsq = h1*h1 + h2*h2 + h3*h3;
        term = -pterm * hsq;
        expterm = 0.0;
        if (term > -50.0 && hsq != 0.0) {
          denom = volterm*hsq*bsmod1[k1]*bsmod2[k2]*bsmod3[k3];
          expterm = exp(term) / denom;
          struc2 = qgrid[k3][k2][k1][0]*qgrip[k3][k2][k1][0] +
            qgrid[k3][k2][k1][1]*qgrip[k3][k2][k1][1];
          eterm = 0.5 * felec * expterm * struc2;
          vterm = (2.0/hsq) * (1.0-term) * eterm;
          virpolar[0] -= h1*h1*vterm - eterm;
          virpolar[1] -= h2*h2*vterm - eterm;
          virpolar[2] -= h3*h3*vterm - eterm;
          virpolar[3] -= h1*h2*vterm;
          virpolar[4] -= h1*h3*vterm;
          virpolar[5] -= h2*h3*vterm;
        }
      }

      // now do the TCG terms with "UBD" dotted with "UAP"

      for (i = 0; i < nlocal; i++) {
        for (j = 0; j < 10; j++)
          cmp[i][j] = 0.0;
        for (j = 1; j < 4; j++)
          cmp[i][j] = uap[i][m][j-1];
      }

      cmp_to_fmp(cmp,fmp);
      grid_mpole(fmp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      for (k = 0; k < nfft3; k++) {
        for (j = 0; j < nfft2; j++) {
          for (i = 0; i < nfft1; i++) {
            qgrip[k][j][i][0] = qgrid[k][j][i][0];
            qgrip[k][j][i][1] = qgrid[k][j][i][1];
          }
        }
      }

      for (i = 0; i < nlocal; i++) {
        for (j = 1; j < 4; j++)
          cmp[i][j] = ubd[i][m][j-1];
      }

      cmp_to_fmp(cmp,fmp);
      grid_mpole(fmp);
      efft->compute(qgrid[0][0][0],qgrid[0][0][0],1);

      // make the scalar summation over reciprocal lattice
      // NOTE: this loop has to be distributed for parallel
      // NOTE: why does this one include m = 0 ?

      for (m = 1; m < ntot; m++) {
        k1 = m % nfft1;
        k2 = (m % nff) / nfft1;
        k3 = m/nff;
        r1 = (k1 >= nf1) ? k1-nfft1 : k1;
        r2 = (k2 >= nf2) ? k2-nfft2 : k2;
        r3 = (k3 >= nf3) ? k3-nfft3 : k3;
        h1 = recip[0][0]*r1 + recip[0][1]*r2 + recip[0][2]*r3;
        h2 = recip[1][0]*r1 + recip[1][1]*r2 + recip[1][2]*r3;
        h3 = recip[2][0]*r1 + recip[2][1]*r2 + recip[2][2]*r3;
        hsq = h1*h1 + h2*h2 + h3*h3;
        term = -pterm * hsq;
        expterm = 0.0;
        if (term > -50.0 && hsq != 0.0) {
          denom = volterm*hsq*bsmod1[k1]*bsmod2[k2]*bsmod3[k3];
          expterm = exp(term) / denom;
          struc2 = qgrid[k3][k2][k1][0]*qgrip[k3][k2][k1][0] +
            qgrid[k3][k2][k1][1]*qgrip[k3][k2][k1][1];
          eterm = 0.5 * felec * expterm * struc2;
          vterm = (2.0/hsq) * (1.0-term) * eterm;
          virpolar[0] -= h1*h1*vterm - eterm;
          virpolar[1] -= h2*h2*vterm - eterm;
          virpolar[2] -= h3*h3*vterm - eterm;
          virpolar[3] -= h1*h2*vterm;
          virpolar[4] -= h1*h3*vterm;
          virpolar[5] -= h2*h3*vterm;
        }
      }
    }
  }
  */

  // increment the total internal virial tensor components

  if (vflag_global) {
    virpolar[0] -= vxx;
    virpolar[1] -= vyy;
    virpolar[2] -= vzz;
    virpolar[3] -= vxy;
    virpolar[4] -= vxz;
    virpolar[5] -= vyz;
  }
}
