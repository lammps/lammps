// **************************************************************************
//                                   amoeba.cu
//                             -------------------
//                          Trung Dac Nguyen (Northwestern)
//
//  Device code for acceleration of the amoeba pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : trung.nguyen@northwestern.edu
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include <stdio.h>
#include "lal_aux_fun1.h"
#ifdef LAMMPS_SMALLBIG
#define tagint int
#endif
#ifdef LAMMPS_BIGBIG
#include "inttypes.h"
#define tagint int64_t
#endif
#ifdef LAMMPS_SMALLSMALL
#define tagint int
#endif
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( q_tex,float);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
#endif

#else
#define pos_tex x_
#define q_tex q_
#endif

#if (SHUFFLE_AVAIL == 0)

#define local_allocate_store_ufld()                                         \
    __local acctyp red_acc[6][BLOCK_PAIR];

#define store_answers_tep(ufld, dufld, ii, inum,tid, t_per_atom, offset,    \
                          i, tep)                                           \
  if (t_per_atom>1) {                                                       \
    red_acc[0][tid]=ufld[0];                                                \
    red_acc[1][tid]=ufld[1];                                                \
    red_acc[2][tid]=ufld[2];                                                \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<3; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    ufld[0]=red_acc[0][tid];                                                \
    ufld[1]=red_acc[1][tid];                                                \
    ufld[2]=red_acc[2][tid];                                                \
    red_acc[0][tid]=dufld[0];                                               \
    red_acc[1][tid]=dufld[1];                                               \
    red_acc[2][tid]=dufld[2];                                               \
    red_acc[3][tid]=dufld[3];                                               \
    red_acc[4][tid]=dufld[4];                                               \
    red_acc[5][tid]=dufld[5];                                               \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    dufld[0]=red_acc[0][tid];                                               \
    dufld[1]=red_acc[1][tid];                                               \
    dufld[2]=red_acc[2][tid];                                               \
    dufld[3]=red_acc[3][tid];                                               \
    dufld[4]=red_acc[4][tid];                                               \
    dufld[5]=red_acc[5][tid];                                               \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    numtyp4 t;                                                              \
    t.x = diz*ufld[1] - diy*ufld[2] + qixz*dufld[1] - qixy*dufld[3] +       \
      (numtyp)2.0*qiyz*(dufld[2]-dufld[5]) + (qizz-qiyy)*dufld[4];          \
    t.y = dix*ufld[2] - diz*ufld[0] - qiyz*dufld[1] + qixy*dufld[4] +       \
      (numtyp)2.0*qixz*(dufld[5]-dufld[0]) + (qixx-qizz)*dufld[3];          \
    t.z = diy*ufld[0] - dix*ufld[1] + qiyz*dufld[3] - qixz*dufld[4] +       \
      (numtyp)2.0*qixy*(dufld[0]-dufld[2]) + (qiyy-qixx)*dufld[1];          \
    tep[i]=t;                                                               \
  }

#define store_answers_fieldp(_fieldp, ii, inum,tid, t_per_atom, offset, i,  \
                              fieldp)                                       \
  if (t_per_atom>1) {                                                       \
    red_acc[0][tid]=_fieldp[0];                                             \
    red_acc[1][tid]=_fieldp[1];                                             \
    red_acc[2][tid]=_fieldp[2];                                             \
    red_acc[3][tid]=_fieldp[3];                                             \
    red_acc[4][tid]=_fieldp[4];                                             \
    red_acc[5][tid]=_fieldp[5];                                             \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    _fieldp[0]=red_acc[0][tid];                                             \
    _fieldp[1]=red_acc[1][tid];                                             \
    _fieldp[2]=red_acc[2][tid];                                             \
    _fieldp[3]=red_acc[3][tid];                                             \
    _fieldp[4]=red_acc[4][tid];                                             \
    _fieldp[5]=red_acc[5][tid];                                             \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    numtyp4 f, fp;                                                          \
    f.x  = _fieldp[0]; f.y  = _fieldp[0]; f.z  = _fieldp[2];                \
    fp.x = _fieldp[3]; fp.y = _fieldp[4]; fp.z = _fieldp[5];                \
    fieldp[ii] = f;                                                         \
    fieldp[ii+inum] = fp;                                                   \
  }

#else

#define local_allocate_store_ufld()

#define store_answers_tep(ufld, dufld, ii, inum,tid, t_per_atom, offset,    \
                          i, tep)                                           \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      ufld[0] += shfl_down(ufld[0], s, t_per_atom);                         \
      ufld[1] += shfl_down(ufld[1], s, t_per_atom);                         \
      ufld[2] += shfl_down(ufld[2], s, t_per_atom);                         \
      dufld[0] += shfl_down(dufld[0], s, t_per_atom);                       \
      dufld[1] += shfl_down(dufld[1], s, t_per_atom);                       \
      dufld[2] += shfl_down(dufld[2], s, t_per_atom);                       \
      dufld[3] += shfl_down(dufld[3], s, t_per_atom);                       \
      dufld[4] += shfl_down(dufld[4], s, t_per_atom);                       \
      dufld[5] += shfl_down(dufld[5], s, t_per_atom);                       \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    numtyp4 t;                                                              \
    t.x = diz*ufld[1] - diy*ufld[2] + qixz*dufld[1] - qixy*dufld[3] +       \
      (numtyp)2.0*qiyz*(dufld[2]-dufld[5]) + (qizz-qiyy)*dufld[4];          \
    t.y = dix*ufld[2] - diz*ufld[0] - qiyz*dufld[1] + qixy*dufld[4] +       \
      (numtyp)2.0*qixz*(dufld[5]-dufld[0]) + (qixx-qizz)*dufld[3];          \
    t.z = diy*ufld[0] - dix*ufld[1] + qiyz*dufld[3] - qixz*dufld[4] +       \
      (numtyp)2.0*qixy*(dufld[0]-dufld[2]) + (qiyy-qixx)*dufld[1];          \
    tep[i]=t;                                                               \
  }

#define store_answers_fieldp(_fieldp, ii, inum,tid, t_per_atom, offset, i,  \
                             fieldp)                                        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      _fieldp[0] += shfl_down(_fieldp[0], s, t_per_atom);                   \
      _fieldp[1] += shfl_down(_fieldp[1], s, t_per_atom);                   \
      _fieldp[2] += shfl_down(_fieldp[2], s, t_per_atom);                   \
      _fieldp[3] += shfl_down(_fieldp[3], s, t_per_atom);                   \
      _fieldp[4] += shfl_down(_fieldp[4], s, t_per_atom);                   \
      _fieldp[5] += shfl_down(_fieldp[5], s, t_per_atom);                   \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    numtyp4 f, fp;                                                          \
    f.x  = _fieldp[0]; f.y  = _fieldp[0]; f.z  = _fieldp[2];                \
    fp.x = _fieldp[3]; fp.y = _fieldp[4]; fp.z = _fieldp[5];                \
    fieldp[ii] = f;                                                         \
    fieldp[ii+inum] = fp;                                                   \
  }

#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MY_PIS (acctyp)1.77245385090551602729

/* ----------------------------------------------------------------------
   polar_real = real-space portion of induced dipole polarization
   adapted from Tinker epreal1d() routine
------------------------------------------------------------------------- */

__kernel void k_amoeba_polar(const __global numtyp4 *restrict x_,
                            const __global numtyp *restrict extra,
                            const __global numtyp4 *restrict damping,
                            const __global numtyp4 *restrict sp_polar,
                            const __global int *dev_nbor,
                            const __global int *dev_packed,
                            __global acctyp4 *restrict ans,
                            __global acctyp *restrict engv,
                            __global numtyp4 *restrict tep,
                            const int eflag, const int vflag, const int inum,
                            const int nall, const int nbor_pitch, const int t_per_atom,
                            const numtyp aewald, const numtyp felec,
                            const numtyp off2, const numtyp polar_dscale,
                            const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_ufld();
  local_allocate_store_charge();

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int l=0; l<6; l++) virial[l]=(acctyp)0;
  }

  acctyp ufld[3];
  ufld[0] = (acctyp)0; ufld[1]=(acctyp)0; ufld[2]=(acctyp)0;
  acctyp dufld[6];
  for (int l=0; l<6; l++) dufld[l]=(acctyp)0;

  numtyp dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz;
  numtyp4* polar1 = (numtyp4*)(&extra[0]);
  numtyp4* polar2 = (numtyp4*)(&extra[4*nall]);
  numtyp4* polar3 = (numtyp4*)(&extra[8*nall]);
  numtyp4* polar4 = (numtyp4*)(&extra[12*nall]);
  numtyp4* polar5 = (numtyp4*)(&extra[16*nall]);

  //numtyp4 xi__;

  if (ii<inum) {
    int k,m,itype,igroup;
    numtyp bfac;
    numtyp psc3,psc5,psc7;
    numtyp dsc3,dsc5,dsc7;
    numtyp usc3,usc5;
    numtyp psr3,psr5,psr7;
    numtyp dsr3,dsr5,dsr7;
    numtyp usr5;
    numtyp term1,term2,term3;
    numtyp term4,term5;
    numtyp term6,term7;
    numtyp rc3[3],rc5[3],rc7[3];
    numtyp prc3[3],prc5[3],prc7[3];
    numtyp drc3[3],drc5[3],drc7[3];
    numtyp urc3[3],urc5[3];
    numtyp bn[5];
    numtyp ci,uix,uiy,uiz,uixp,uiyp,uizp;

    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    //numtyp qtmp; fetch(qtmp,i,q_tex);
    //int itype=ix.w;

    ci  = polar1[i].x;    // rpole[i][0];
    dix = polar1[i].y;    // rpole[i][1];
    diy = polar1[i].z;    // rpole[i][2];
    diz = polar1[i].w;    // rpole[i][3];
    qixx = polar2[i].x;   // rpole[i][4];
    qixy = polar2[i].y;   // rpole[i][5];
    qixz = polar2[i].z;   // rpole[i][6];
    qiyy = polar2[i].w;   // rpole[i][8];
    qiyz   = polar3[i].x; // rpole[i][9];
    qizz   = polar3[i].y; // rpole[i][12];
    itype  = polar3[i].z; // amtype[i];
    igroup = polar3[i].w; // amgroup[i];
    uix = polar4[i].x;    // uind[i][0];
    uiy = polar4[i].y;    // uind[i][1];
    uiz = polar4[i].z;    // uind[i][2];
    uixp = polar5[i].x;   // uinp[i][0];
    uiyp = polar5[i].y;   // uinp[i][1];
    uizp = polar5[i].z;   // uinp[i][2];

    // debug:
    // xi__ = ix; xi__.w = itype;

    numtyp pdi = damping[itype].x;
    numtyp pti = damping[itype].y;

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=dev_packed[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      //int jtype=jx.w;
 
      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      if (r2>off2) continue;
  
      numtyp r = ucl_sqrt(r2);
      
      numtyp ck = polar1[j].x;   // rpole[j][0];
      numtyp dkx = polar1[j].y;  // rpole[j][1];
      numtyp dky = polar1[j].z;  // rpole[j][2];
      numtyp dkz = polar1[j].w;  // rpole[j][3];
      numtyp qkxx = polar2[j].x; // rpole[j][4];
      numtyp qkxy = polar2[j].y; // rpole[j][5];
      numtyp qkxz = polar2[j].z; // rpole[j][6];
      numtyp qkyy = polar2[j].w; // rpole[j][8];
      numtyp qkyz = polar3[j].x; // rpole[j][9];
      numtyp qkzz = polar3[j].y; // rpole[j][12];
      int jtype =   polar3[j].z; // amtype[j];
      int jgroup =  polar3[j].w; // amgroup[j];
      numtyp ukx = polar4[j].x;  // uind[j][0];
      numtyp uky = polar4[j].y;  // uind[j][1];
      numtyp ukz = polar4[j].z;  // uind[j][2];
      numtyp ukxp = polar5[j].x; // uinp[j][0];
      numtyp ukyp = polar5[j].y; // uinp[j][1];
      numtyp ukzp = polar5[j].z; // uinp[j][2];

      numtyp factor_wscale, factor_dscale, factor_pscale, factor_uscale;
      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      factor_wscale = sp_pol.x; // sp_polar_wscale[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_pscale = sp_pol.y; // sp_polar_piscale[sbmask15(jextra)];
        factor_dscale = polar_dscale;
        factor_uscale = polar_uscale;
      } else {
        factor_pscale = sp_pol.z; // sp_polar_pscale[sbmask15(jextra)];
        factor_dscale = factor_uscale = (numtyp)1.0;
      }

      // intermediates involving moments and separation distance

      numtyp dir = dix*xr + diy*yr + diz*zr;
      numtyp qix = qixx*xr + qixy*yr + qixz*zr;
      numtyp qiy = qixy*xr + qiyy*yr + qiyz*zr;
      numtyp qiz = qixz*xr + qiyz*yr + qizz*zr;
      numtyp qir = qix*xr + qiy*yr + qiz*zr;
      numtyp dkr = dkx*xr + dky*yr + dkz*zr;
      numtyp qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      numtyp qky = qkxy*xr + qkyy*yr + qkyz*zr;
      numtyp qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      numtyp qkr = qkx*xr + qky*yr + qkz*zr;
      numtyp uir = uix*xr + uiy*yr + uiz*zr;
      numtyp uirp = uixp*xr + uiyp*yr + uizp*zr;
      numtyp ukr = ukx*xr + uky*yr + ukz*zr;
      numtyp ukrp = ukxp*xr + ukyp*yr + ukzp*zr;

      // get reciprocal distance terms for this interaction

      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = felec * rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;
      numtyp rr9 = (numtyp)7.0 * rr7 * r2inv;

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a;
      //bn[0] = erfc(ralpha) / r;
      bn[0] = _erfc * rinv;
      numtyp alsq2 = (numtyp)2.0 * aewald*aewald;
      numtyp alsq2n = (numtyp)0.0;
      if (aewald > (numtyp)0.0) alsq2n = (numtyp)1.0 / (MY_PIS*aewald);

      for (m = 1; m <= 4; m++) {
        bfac = (numtyp) (m+m-1);
        alsq2n = alsq2 * alsq2n;
        bn[m] = (bfac*bn[m-1]+alsq2n*exp2a) / r2;
      }
      for (m = 0; m < 5; m++) bn[m] *= felec;

      // apply Thole polarization damping to scale factors

      numtyp sc3 = (numtyp)1.0;
      numtyp sc5 = (numtyp)1.0;
      numtyp sc7 = (numtyp)1.0;
      for (k = 0; k < 3; k++) {
        rc3[k] = (numtyp)0.0;
        rc5[k] = (numtyp)0.0;
        rc7[k] = (numtyp)0.0;
      }

      // apply Thole polarization damping to scale factors

      numtyp damp = pdi * damping[jtype].x; // pdamp[jtype]
      if (damp != (numtyp)0.0) {
        numtyp pgamma = MIN(pti,damping[jtype].y); // thole[jtype]
        damp = pgamma * ucl_powr(r/damp,(numtyp)3.0);
        if (damp < (numtyp)50.0) {
          numtyp expdamp = ucl_exp(-damp);
          sc3 = (numtyp)1.0 - expdamp;
          sc5 = (numtyp)1.0 - ((numtyp)1.0+damp)*expdamp;
          sc7 = (numtyp)1.0 - ((numtyp)1.0+damp+(numtyp)0.6*damp*damp) * expdamp;
          numtyp temp3 = (numtyp)3.0 * damp * expdamp * r2inv;
          numtyp temp5 = damp;
          numtyp temp7 = (numtyp)-0.2 + (numtyp)0.6*damp;
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

        psc3 = (numtyp)1.0 - sc3*factor_pscale;
        psc5 = (numtyp)1.0 - sc5*factor_pscale;
        psc7 = (numtyp)1.0 - sc7*factor_pscale;
        dsc3 = (numtyp)1.0 - sc3*factor_dscale;
        dsc5 = (numtyp)1.0 - sc5*factor_dscale;
        dsc7 = (numtyp)1.0 - sc7*factor_dscale;
        usc3 = (numtyp)1.0 - sc3*factor_uscale;
        usc5 = (numtyp)1.0 - sc5*factor_uscale;
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
      } else { // damp == 0: ???
      }

      // get the induced dipole field used for dipole torques

      numtyp tix3 = psr3*ukx + dsr3*ukxp;
      numtyp tiy3 = psr3*uky + dsr3*ukyp;
      numtyp tiz3 = psr3*ukz + dsr3*ukzp;
      numtyp tuir = -psr5*ukr - dsr5*ukrp;
      
      ufld[0] += tix3 + xr*tuir;
      ufld[1] += tiy3 + yr*tuir;
      ufld[2] += tiz3 + zr*tuir;

      // get induced dipole field gradient used for quadrupole torques

      numtyp tix5 = (numtyp)2.0 * (psr5*ukx+dsr5*ukxp);
      numtyp tiy5 = (numtyp)2.0 * (psr5*uky+dsr5*ukyp);
      numtyp tiz5 = (numtyp)2.0 * (psr5*ukz+dsr5*ukzp);
      tuir = -psr7*ukr - dsr7*ukrp;
      
      dufld[0] += xr*tix5 + xr*xr*tuir;
      dufld[1] += xr*tiy5 + yr*tix5 + (numtyp)2.0*xr*yr*tuir;
      dufld[2] += yr*tiy5 + yr*yr*tuir;
      dufld[3] += xr*tiz5 + zr*tix5 + (numtyp)2.0*xr*zr*tuir;
      dufld[4] += yr*tiz5 + zr*tiy5 + (numtyp)2.0*yr*zr*tuir;
      dufld[5] += zr*tiz5 + zr*zr*tuir;
      
      // get the dEd/dR terms used for direct polarization force

      term1 = bn[2] - dsc3*rr5;
      term2 = bn[3] - dsc5*rr7;
      term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3[0];
      term4 = rr3*drc3[0] - term1*xr - dsr5*xr;
      term5 = term2*xr*xr - dsr5 - rr5*xr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*xr - bn[3] - rr7*xr*drc7[0];
      term7 = rr5*drc5[0] - (numtyp)2.0*bn[3]*xr + (dsc5+(numtyp)1.5*dsc7)*rr7*xr;
      numtyp tixx = ci*term3 + dix*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7 + (numtyp)2.0*qix*term7 + qir*term6;
      numtyp tkxx = ck*term3 - dkx*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7 + (numtyp)2.0*qkx*term7 + qkr*term6;

      term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3[1];
      term4 = rr3*drc3[1] - term1*yr - dsr5*yr;
      term5 = term2*yr*yr - dsr5 - rr5*yr*drc5[1];
      term6 = (bn[4]-dsc7*rr9)*yr*yr - bn[3] - rr7*yr*drc7[1];
      term7 = rr5*drc5[1] - (numtyp)2.0*bn[3]*yr + (dsc5+(numtyp)1.5*dsc7)*rr7*yr;
      numtyp tiyy = ci*term3 + diy*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7 + (numtyp)2.0*qiy*term7 + qir*term6;
      numtyp tkyy = ck*term3 - dky*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7 + (numtyp)2.0*qky*term7 + qkr*term6;

      term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3[2];
      term4 = rr3*drc3[2] - term1*zr - dsr5*zr;
      term5 = term2*zr*zr - dsr5 - rr5*zr*drc5[2];
      term6 = (bn[4]-dsc7*rr9)*zr*zr - bn[3] - rr7*zr*drc7[2];
      term7 = rr5*drc5[2] - (numtyp)2.0*bn[3]*zr + (dsc5+(numtyp)1.5*dsc7)*rr7*zr;
      numtyp tizz = ci*term3 + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7 + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkzz = ck*term3 - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7 + (numtyp)2.0*qkz*term7 + qkr*term6;

      term3 = term1*xr*yr - rr3*yr*drc3[0];
      term4 = rr3*drc3[0] - term1*xr;
      term5 = term2*xr*yr - rr5*yr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*yr - rr7*yr*drc7[0];
      term7 = rr5*drc5[0] - term2*xr;
      numtyp tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixy - (numtyp)2.0*dsr7*yr*qix + (numtyp)2.0*qiy*term7 + qir*term6;
      numtyp tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxy - (numtyp)2.0*dsr7*yr*qkx +(numtyp) 2.0*qky*term7 + qkr*term6;

      term3 = term1*xr*zr - rr3*zr*drc3[0];
      term5 = term2*xr*zr - rr5*zr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*zr - rr7*zr*drc7[0];
      numtyp tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixz - (numtyp)2.0*dsr7*zr*qix + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxz - (numtyp)2.0*dsr7*zr*qkx + (numtyp)2.0*qkz*term7 + qkr*term6;

      term3 = term1*yr*zr - rr3*zr*drc3[1];
      term4 = rr3*drc3[1] - term1*yr;
      term5 = term2*yr*zr - rr5*zr*drc5[1];
      term6 = (bn[4]-dsc7*rr9)*yr*zr - rr7*zr*drc7[1];
      term7 = rr5*drc5[1] - term2*yr;
      numtyp tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qiyz - (numtyp)2.0*dsr7*zr*qiy + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkyz - (numtyp)2.0*dsr7*zr*qky + (numtyp)2.0*qkz*term7 + qkr*term6;

      numtyp depx = tixx*ukxp + tixy*ukyp + tixz*ukzp - tkxx*uixp - tkxy*uiyp - tkxz*uizp;
      numtyp depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp - tkxy*uixp - tkyy*uiyp - tkyz*uizp;
      numtyp depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp - tkxz*uixp - tkyz*uiyp - tkzz*uizp;

      numtyp frcx = depx;
      numtyp frcy = depy;
      numtyp frcz = depz;

      // get the dEp/dR terms used for direct polarization force
      
      // tixx and tkxx
      term1 = bn[2] - psc3*rr5;
      term2 = bn[3] - psc5*rr7;
      term3 = -psr3 + term1*xr*xr - rr3*xr*prc3[0];
      term4 = rr3*prc3[0] - term1*xr - psr5*xr;
      term5 = term2*xr*xr - psr5 - rr5*xr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*xr - bn[3] - rr7*xr*prc7[0];
      term7 = rr5*prc5[0] - (numtyp)2.0*bn[3]*xr + (psc5+(numtyp)1.5*psc7)*rr7*xr;
      tixx = ci*term3 + dix*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7 + (numtyp)2.0*qix*term7 + qir*term6;
      tkxx = ck*term3 - dkx*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7 + (numtyp)2.0*qkx*term7 + qkr*term6;

      // tiyy and tkyy
      term3 = -psr3 + term1*yr*yr - rr3*yr*prc3[1];
      term4 = rr3*prc3[1] - term1*yr - psr5*yr;
      term5 = term2*yr*yr - psr5 - rr5*yr*prc5[1];
      term6 = (bn[4]-psc7*rr9)*yr*yr - bn[3] - rr7*yr*prc7[1];
      term7 = rr5*prc5[1] - (numtyp)2.0*bn[3]*yr + (psc5+(numtyp)1.5*psc7)*rr7*yr;
      tiyy = ci*term3 + diy*term4 + dir*term5 +
        (numtyp)2.0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7 + (numtyp)2.0*qiy*term7 + qir*term6;
      tkyy = ck*term3 - dky*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7 + (numtyp)2.0*qky*term7 + qkr*term6;

      // tizz and tkzz
      term3 = -psr3 + term1*zr*zr - rr3*zr*prc3[2];
      term4 = rr3*prc3[2] - term1*zr - psr5*zr;
      term5 = term2*zr*zr - psr5 - rr5*zr*prc5[2];
      term6 = (bn[4]-psc7*rr9)*zr*zr - bn[3] - rr7*zr*prc7[2];
      term7 = rr5*prc5[2] - (numtyp)2.0*bn[3]*zr + (psc5+(numtyp)1.5*psc7)*rr7*zr;
      tizz = ci*term3 + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7 + (numtyp)2.0*qiz*term7 + qir*term6;
      tkzz = ck*term3 - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7 + (numtyp)2.0*qkz*term7 + qkr*term6;

      // tixy and tkxy
      term3 = term1*xr*yr - rr3*yr*prc3[0];
      term4 = rr3*prc3[0] - term1*xr;
      term5 = term2*xr*yr - rr5*yr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*yr - rr7*yr*prc7[0];
      term7 = rr5*prc5[0] - term2*xr;
      tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixy - (numtyp)2.0*psr7*yr*qix + (numtyp)2.0*qiy*term7 + qir*term6;
      tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxy - (numtyp)2.0*psr7*yr*qkx + (numtyp)2.0*qky*term7 + qkr*term6;

      // tixz and tkxz
      term3 = term1*xr*zr - rr3*zr*prc3[0];
      term5 = term2*xr*zr - rr5*zr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*zr - rr7*zr*prc7[0];
      tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixz - (numtyp)2.0*psr7*zr*qix + (numtyp)2.0*qiz*term7 + qir*term6;
      tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxz - (numtyp)2.0*psr7*zr*qkx + (numtyp)2.0*qkz*term7 + qkr*term6;

      // tiyz and tkyz
      term3 = term1*yr*zr - rr3*zr*prc3[1];
      term4 = rr3*prc3[1] - term1*yr;
      term5 = term2*yr*zr - rr5*zr*prc5[1];
      term6 = (bn[4]-psc7*rr9)*yr*zr - rr7*zr*prc7[1];
      term7 = rr5*prc5[1] - term2*yr;
      tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qiyz - (numtyp)2.0*psr7*zr*qiy + (numtyp)2.0*qiz*term7 + qir*term6;
      tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkyz - (numtyp)2.0*psr7*zr*qky + (numtyp)2.0*qkz*term7 + qkr*term6;

      depx = tixx*ukx + tixy*uky + tixz*ukz - tkxx*uix - tkxy*uiy - tkxz*uiz;
      depy = tixy*ukx + tiyy*uky + tiyz*ukz - tkxy*uix - tkyy*uiy - tkyz*uiz;
      depz = tixz*ukx + tiyz*uky + tizz*ukz - tkxz*uix - tkyz*uiy - tkzz*uiz;

      frcx = frcx + depx;
      frcy = frcy + depy;
      frcz = frcz + depz;

      // get the dtau/dr terms used for mutual polarization force
      // poltyp == MUTUAL  && amoeba
          
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

      f.x -= frcx;
      f.y -= frcy;
      f.z -= frcz;

      if (EVFLAG && vflag) {
        numtyp vxx = xr * frcx;
        numtyp vxy = (numtyp)0.5 * (yr*frcx+xr*frcy);
        numtyp vxz = (numtyp)0.5 * (zr*frcx+xr*frcz);
        numtyp vyy = yr * frcy;
        numtyp vyz = (numtyp)0.5 * (zr*frcy+yr*frcz);
        numtyp vzz = zr * frcz;

        virial[0] += vxx;
        virial[1] += vyy;
        virial[2] += vzz;
        virial[3] += vxy;
        virial[4] += vxz;
        virial[5] += vyz;
      }
    } // nbor
    
  } // ii<inum

  // accumulate ufld and dufld to compute tep
  store_answers_tep(ufld,dufld,ii,inum,tid,t_per_atom,offset,i,tep);

  // accumate force, energy and virial
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv);
}

/* ----------------------------------------------------------------------
  udirect2b = Ewald real direct field via list
  udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

__kernel void k_amoeba_udirect2b(const __global numtyp4 *restrict x_,
                                 const __global numtyp *restrict extra,
                                 const __global numtyp4 *restrict damping,
                                 const __global numtyp4 *restrict sp_polar,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 __global numtyp4 *restrict fieldp,
                                 const int inum,  const int nall,
                                 const int nbor_pitch, const int t_per_atom,
                                 const numtyp aewald, const numtyp off2,
                                 const numtyp polar_dscale, const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_charge();

  acctyp _fieldp[6];
  for (int l=0; l<6; l++) _fieldp[l]=(acctyp)0;

  numtyp dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz;
  numtyp4* polar1 = (numtyp4*)(&extra[0]);
  numtyp4* polar2 = (numtyp4*)(&extra[4*nall]);
  numtyp4* polar3 = (numtyp4*)(&extra[8*nall]);

  //numtyp4 xi__;

  if (ii<inum) {
    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    //numtyp qtmp; fetch(qtmp,i,q_tex);
    //int itype=ix.w;

    int itype,igroup;
    numtyp bn[4],bcn[3];
    numtyp fid[3],fip[3];
    numtyp ci,uix,uiy,uiz,uixp,uiyp,uizp;
    
    ci  = polar1[i].x;    // rpole[i][0];
    dix = polar1[i].y;    // rpole[i][1];
    diy = polar1[i].z;    // rpole[i][2];
    diz = polar1[i].w;    // rpole[i][3];
    qixx = polar2[i].x;   // rpole[i][4];
    qixy = polar2[i].y;   // rpole[i][5];
    qixz = polar2[i].z;   // rpole[i][6];
    qiyy = polar2[i].w;   // rpole[i][8];
    qiyz   = polar3[i].x; // rpole[i][9];
    qizz   = polar3[i].y; // rpole[i][12];
    itype  = polar3[i].z; // amtype[i];
    igroup = polar3[i].w; // amgroup[i];
    
    // debug:
    // xi__ = ix; xi__.w = itype;

    numtyp pdi = damping[itype].x;
    numtyp pti = damping[itype].y;
    numtyp ddi = damping[itype].z;

    numtyp aesq2 = (numtyp)2.0 * aewald*aewald;
    numtyp aesq2n = (numtyp)0.0;
    if (aewald > (numtyp)0.0) aesq2n = (numtyp)1.0 / (MY_PIS*aewald);

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=dev_packed[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      //int jtype=jx.w;
 
      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      if (r2>off2) continue;
  
      numtyp r = ucl_sqrt(r2);
      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;

      numtyp ck = polar1[j].x;   // rpole[j][0];
      numtyp dkx = polar1[j].y;  // rpole[j][1];
      numtyp dky = polar1[j].z;  // rpole[j][2];
      numtyp dkz = polar1[j].w;  // rpole[j][3];
      numtyp qkxx = polar2[j].x; // rpole[j][4];
      numtyp qkxy = polar2[j].y; // rpole[j][5];
      numtyp qkxz = polar2[j].z; // rpole[j][6];
      numtyp qkyy = polar2[j].w; // rpole[j][8];
      numtyp qkyz = polar3[j].x; // rpole[j][9];
      numtyp qkzz = polar3[j].y; // rpole[j][12];
      int jtype =   polar3[j].z; // amtype[j];
      int jgroup =  polar3[j].w; // amgroup[j];

      numtyp factor_wscale, factor_dscale, factor_pscale, factor_uscale;
      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      factor_wscale = sp_pol.x; // sp_polar_wscale[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_pscale = sp_pol.y; // sp_polar_piscale[sbmask15(jextra)];
        factor_dscale = polar_dscale;
        factor_uscale = polar_uscale;
      } else {
        factor_pscale = sp_pol.z; // sp_polar_pscale[sbmask15(jextra)];
        factor_dscale = factor_uscale = (numtyp)1.0;
      }

      // intermediates involving moments and separation distance

      numtyp dir = dix*xr + diy*yr + diz*zr;
      numtyp qix = qixx*xr + qixy*yr + qixz*zr;
      numtyp qiy = qixy*xr + qiyy*yr + qiyz*zr;
      numtyp qiz = qixz*xr + qiyz*yr + qizz*zr;
      numtyp qir = qix*xr + qiy*yr + qiz*zr;
      numtyp dkr = dkx*xr + dky*yr + dkz*zr;
      numtyp qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      numtyp qky = qkxy*xr + qkyy*yr + qkyz*zr;
      numtyp qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      numtyp qkr = qkx*xr + qky*yr + qkz*zr;

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a;
      //bn[0] = erfc(ralpha) / r;
      bn[0] = _erfc * rinv;

      numtyp aefac = aesq2n;
      for (int m = 1; m <= 3; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * r2inv;
      }

      // find the field components for Thole polarization damping

      numtyp scale3 = (numtyp)1.0;
      numtyp scale5 = (numtyp)1.0;
      numtyp scale7 = (numtyp)1.0;
      numtyp damp = pdi * damping[jtype].x; // pdamp[jtype]
      if (damp != (numtyp)0.0) {
        numtyp pgamma = MIN(ddi,damping[jtype].z); // dirdamp[jtype]
        if (pgamma != (numtyp)0.0) {
          damp = pgamma * ucl_powr(r/damp,(numtyp)1.5);
          if (damp < (numtyp)50.0) {
            numtyp expdamp = ucl_exp(-damp) ;
            scale3 = (numtyp)1.0 - expdamp ;
            scale5 = (numtyp)1.0 - expdamp*((numtyp)1.0+(numtyp)0.5*damp);
            scale7 = (numtyp)1.0 - expdamp*((numtyp)1.0+(numtyp)0.65*damp + (numtyp)0.15*damp*damp);
          }
        } else {
          pgamma = MIN(pti,damping[jtype].y); // thole[jtype]
          damp = pgamma * ucl_powr(r/damp,3.0);
          if (damp < (numtyp)50.0) {
            numtyp expdamp = ucl_exp(-damp);
            scale3 = (numtyp)1.0 - expdamp;
            scale5 = (numtyp)1.0 - expdamp*((numtyp)1.0+damp);
            scale7 = (numtyp)1.0 - expdamp*((numtyp)1.0+damp + (numtyp)0.6*damp*damp);
          }
        }
      } else { // damp == 0: ???
      }

      numtyp scalek = factor_dscale;
      bcn[0] = bn[1] - ((numtyp)1.0-scalek*scale3)*rr3;
      bcn[1] = bn[2] - ((numtyp)1.0-scalek*scale5)*rr5;
      bcn[2] = bn[3] - ((numtyp)1.0-scalek*scale7)*rr7;
      fid[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkx + (numtyp)2.0*bcn[1]*qkx;
      fid[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dky + (numtyp)2.0*bcn[1]*qky;
      fid[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkz + (numtyp)2.0*bcn[1]*qkz;
        
      scalek = factor_pscale;
      bcn[0] = bn[1] - ((numtyp)1.0-scalek*scale3)*rr3;
      bcn[1] = bn[2] - ((numtyp)1.0-scalek*scale5)*rr5;
      bcn[2] = bn[3] - ((numtyp)1.0-scalek*scale7)*rr7;
      fip[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkx + (numtyp)2.0*bcn[1]*qkx;
      fip[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dky + (numtyp)2.0*bcn[1]*qky;
      fip[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkz + (numtyp)2.0*bcn[1]*qkz;

      _fieldp[0] += fid[0];
      _fieldp[1] += fid[1];
      _fieldp[2] += fid[2];
      _fieldp[3] += fip[0];
      _fieldp[4] += fip[1];
      _fieldp[5] += fip[2];

    } // nbor

  } // ii<inum

  // accumulate field and fieldp
  
  store_answers_fieldp(_fieldp,ii,inum,tid,t_per_atom,offset,i,fieldp);
}

/* ----------------------------------------------------------------------
   scan standard neighbor list and make it compatible with 1-5 neighbors
   if IJ entry is a 1-2,1-3,1-4 neighbor then adjust offset to SBBITS15
   else scan special15 to see if a 1-5 neighbor and adjust offset to SBBITS15
   else do nothing to IJ entry
------------------------------------------------------------------------- */

__kernel void k_special15(__global int * dev_nbor,
                          const __global int * dev_packed,
                          const __global tagint *restrict tag,
                          const __global int *restrict nspecial15,
                          const __global tagint *restrict special15,
                          const int inum, const int nall, const int nbor_pitch,
                          const int t_per_atom) {
  int tid, ii, offset, n_stride, i;
  atom_info(t_per_atom,ii,tid,offset);

  if (ii<inum) {
  
    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    int n15 = nspecial15[ii];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int sj=dev_packed[nbor];
      int which = sj >> SBBITS & 3;
      int j = sj & NEIGHMASK;
      tagint jtag = tag[j];

      if (!which) {
        int offset=ii;
        for (int k=0; k<n15; k++) {
          if (special15[offset] == jtag) {
            which = 4;
            break;
          }
          offset += nall;
        }
      }

      if (which) dev_nbor[nbor] = j ^ (which << SBBITS15);
    } // for nbor

  } // if ii
}

