// **************************************************************************
//                               lj_tip4p_long.cu
//                             -------------------
//                              V. Nikolskiy (HSE)
//
//  Device code for acceleration of the lj/tip4p/long pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : thevsevak@gmail.com
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)

#include "lal_aux_fun1.h"
#ifdef LAMMPS_SMALLBIG
#define tagint int
#endif
#ifdef LAMMPS_BIGBIG
#ifdef USE_OPENCL
#define tagint long
#else
#include "stdint.h"
#define tagint int64_t
#endif
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
#ifdef LAMMPS_SMALLBIG
#define tagint int
#endif
#ifdef LAMMPS_BIGBIG
#ifdef USE_OPENCL
#define tagint long
#else
#include "stdint.h"
#define tagint int64_t
#endif
#endif
#ifdef LAMMPS_SMALLSMALL
#define tagint int
#endif
#define pos_tex x_
#define q_tex q_
#endif

ucl_inline int atom_mapping(const __global int *map, tagint glob) {
  return map[glob];
}

ucl_inline int closest_image(int i, int j, const __global int* sametag,
                             const __global numtyp4 *restrict x_)
{
  if (j < 0) return j;

  numtyp4 xi; fetch4(xi,i,pos_tex); // = x[i];
  numtyp4 xj; fetch4(xj,j,pos_tex);

  int closest = j;
  numtyp delx = xi.x - xj.x;
  numtyp dely = xi.y - xj.y;
  numtyp delz = xi.z - xj.z;
  numtyp rsqmin = delx*delx + dely*dely + delz*delz;
  numtyp rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    fetch4(xj,j,pos_tex);
    delx = xi.x - xj.x;
    dely = xi.y - xj.y;
    delz = xi.z - xj.z;
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }

  return closest;
}

ucl_inline void compute_newsite(int iO, int  iH1, int  iH2,
    __global numtyp4 *xM, numtyp q,
    numtyp alpha, const __global numtyp4 *restrict x_) {
  numtyp4 xO;  fetch4(xO,iO,pos_tex);
  numtyp4 xH1; fetch4(xH1,iH1,pos_tex);
  numtyp4 xH2; fetch4(xH2,iH2,pos_tex);
  numtyp4 M;

  numtyp delx1 = xH1.x - xO.x;
  numtyp dely1 = xH1.y - xO.y;
  numtyp delz1 = xH1.z - xO.z;

  numtyp delx2 = xH2.x - xO.x;
  numtyp dely2 = xH2.y - xO.y;
  numtyp delz2 = xH2.z - xO.z;

  numtyp ap = alpha * (numtyp)0.5;

  M.x = xO.x + ap * (delx1 + delx2);
  M.y = xO.y + ap * (dely1 + dely2);
  M.z = xO.z + ap * (delz1 + delz2);
  M.w = q;

  *xM = M;
}

__kernel void k_lj_tip4p_long_distrib(const __global numtyp4 *restrict x_,
    __global acctyp4 *restrict ans,
    __global acctyp *restrict engv,
    const int eflag, const int vflag, const int inum,
    const int nbor_pitch, const int t_per_atom,
    __global int *restrict hneigh,
    __global numtyp4 *restrict m,
    const int typeO, const int typeH,
    const numtyp alpha,
    const __global numtyp *restrict q_, const __global acctyp4 *restrict ansO) {

  int i = BLOCK_ID_X*(BLOCK_SIZE_X)+THREAD_ID_X;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;

  if (i<inum) {
    numtyp4 ix; fetch4(ix,i,pos_tex);// = x_[i];
    int itype = ix.w;
    acctyp4 fM, vM;
    acctyp eM;
    // placement of the virial in engv depends on eflag value
    int engv_iter = eflag ? 2 : 0;
    if (itype == typeH) {
      int iO = hneigh[i*4];
      if (iO < inum) {
        fM = ansO[iO];
        f.x += fM.x * (acctyp)0.5 * alpha;
        f.y += fM.y * (acctyp)0.5 * alpha;
        f.z += fM.z * (acctyp)0.5 * alpha;
        if (EVFLAG && vflag) {
          vM = ansO[inum  +iO];
          engv[inum*engv_iter + i] += vM.x * (acctyp)0.5 * alpha; engv_iter++;
          engv[inum*engv_iter + i] += vM.y * (acctyp)0.5 * alpha; engv_iter++;
          engv[inum*engv_iter + i] += vM.z * (acctyp)0.5 * alpha; engv_iter++;
          vM = ansO[inum*2+iO];
          engv[inum*engv_iter + i] += vM.x * (acctyp)0.5 * alpha; engv_iter++;
          engv[inum*engv_iter + i] += vM.y * (acctyp)0.5 * alpha; engv_iter++;
          engv[inum*engv_iter + i] += vM.z * (acctyp)0.5 * alpha;
        }
      }
    }
    if (itype == typeO) {
      fM = ansO[i];
      int iH1 = hneigh[i*4  ];
      int iH2 = hneigh[i*4+1];
      f.x += fM.x * (acctyp)(1 - alpha);
      f.y += fM.y * (acctyp)(1 - alpha);
      f.z += fM.z * (acctyp)(1 - alpha);
      if (EVFLAG && eflag) {
        eM = engv[i+inum];
        engv[inum+i] = eM*(acctyp)(1 - alpha);
        if (iH1 < inum) engv[inum+iH1] += eM * (acctyp)0.5 * alpha;
        if (iH2 < inum) engv[inum+iH2] += eM * (acctyp)0.5 * alpha;
      }
      if (EVFLAG && vflag) {
        vM = ansO[inum   + i];
        engv[inum*engv_iter + i] += vM.x * (acctyp)(1 - alpha); engv_iter++;
        engv[inum*engv_iter + i] += vM.y * (acctyp)(1 - alpha); engv_iter++;
        engv[inum*engv_iter + i] += vM.z * (acctyp)(1 - alpha); engv_iter++;
        vM = ansO[inum*2 + i];
        engv[inum*engv_iter + i] += vM.x * (acctyp)(1 - alpha); engv_iter++;
        engv[inum*engv_iter + i] += vM.y * (acctyp)(1 - alpha); engv_iter++;
        engv[inum*engv_iter + i] += vM.z * (acctyp)(1 - alpha);
      }
    }
    acctyp4 old=ans[i];
    old.x+=f.x;
    old.y+=f.y;
    old.z+=f.z;
    ans[i]=old;
  } // if ii
}

__kernel void k_lj_tip4p_reneigh(const __global numtyp4 *restrict x_,
    const __global int * dev_nbor,
    const __global int * dev_packed,
    const int nall, const int inum,
    const int nbor_pitch, const int t_per_atom,
    __global int *restrict hneigh,
    __global numtyp4 *restrict m,
    const int typeO, const int typeH,
    const __global tagint *restrict tag, const __global int *restrict map,
    const __global int *restrict sametag) {

  int i = BLOCK_ID_X*(BLOCK_SIZE_X)+THREAD_ID_X;

  if (i<nall) {
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    int iH1, iH2, iO;
    int itype = ix.w;
    if(itype == typeO) {
      iO  = i;
      if (hneigh[i*4+2] != -1) {
        iH1 = atom_mapping(map, tag[i] + 1);
        iH2 = atom_mapping(map, tag[i] + 2);
        // set iH1,iH2 to closest image to O
        iH1 = closest_image(i, iH1, sametag, x_);
        iH2 = closest_image(i, iH2, sametag, x_);
        hneigh[i*4  ] = iH1;
        hneigh[i*4+1] = iH2;
        hneigh[i*4+2] = -1;
      }
    }
    if (itype == typeH) {
      if (hneigh[i*4+2] != -1) {
        int iI, iH;
        iI = atom_mapping(map,tag[i] - 1);
        numtyp4 iIx; fetch4(iIx,iI,pos_tex); //x_[iI];
        if ((int)iIx.w == typeH) {
          iO = atom_mapping(map,tag[i] - 2);
          iO  = closest_image(i, iO, sametag, x_);
          iH1 = closest_image(i, iI, sametag, x_);
          iH2 = i;
        } else { //if ((int)iIx.w == typeO)
          iH = atom_mapping(map, tag[i] + 1);
          iO  = closest_image(i,iI,sametag, x_);
          iH1 = i;
          iH2 = closest_image(i,iH,sametag, x_);
        }
        hneigh[i*4+0] = iO;
        hneigh[i*4+1] += -1;
        hneigh[i*4+2] = -1;
      }
    }
  }
}


__kernel void k_lj_tip4p_newsite(const __global numtyp4 *restrict x_,
    const __global int * dev_nbor,
    const __global int * dev_packed,
    const int nall, const int inum,
    const int nbor_pitch, const int t_per_atom,
    __global int *restrict hneigh,
    __global numtyp4 *restrict m,
    const int typeO, const int typeH,
    const numtyp alpha, const __global numtyp *restrict q_) {

  int i = BLOCK_ID_X*(BLOCK_SIZE_X)+THREAD_ID_X;

  if (i<nall) {
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype = ix.w;
    if (itype == typeO) {
      int iH1, iH2, iO;
      iH1 = hneigh[i*4  ];
      iH2 = hneigh[i*4+1];
      iO  = i;
      numtyp qO; fetch(qO,iO,q_tex);
      compute_newsite(iO,iH1,iH2, &m[iO], qO, alpha, x_);
    }
  }
}

__kernel void k_lj_tip4p_long(const __global numtyp4 *restrict x_,
    const __global numtyp4 *restrict lj1,
    const __global numtyp4 *restrict lj3,
    const int lj_types,
    const __global numtyp *restrict sp_lj,
    const __global int * dev_nbor,
    const __global int * dev_packed,
    __global acctyp4 *restrict ans,
    __global acctyp *restrict engv,
    const int eflag, const int vflag, const int inum,
    const int nbor_pitch, const int t_per_atom,
    __global int *restrict hneigh,
    __global numtyp4 *restrict m,
    const int typeO, const int typeH,
    const numtyp alpha,
    const __global numtyp *restrict q_,
    const __global numtyp *restrict cutsq,
    const numtyp qqrd2e, const numtyp g_ewald,
    const numtyp cut_coulsq, const numtyp cut_coulsqplus,
    __global acctyp4 *restrict ansO) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_charge();

  acctyp4 f, fO;
  f.x=(acctyp)0;  f.y=(acctyp)0;  f.z=(acctyp)0;
  fO.x=(acctyp)0; fO.y=(acctyp)0; fO.z=(acctyp)0;
  acctyp energy, e_coul, virial[6], vO[6];
  if (EVFLAG) {
    energy = (acctyp)0;
    e_coul = (acctyp)0;
    for (int i=0; i<6; i++) {
      virial[i]=(acctyp)0;
      vO[i]=(acctyp)0;
    }
  }

  int i;
  if (ii<inum) {
    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype = ix.w;
    numtyp4 x1 = ix;

    int non_local_oxy = 0;
    int iH1, iH2, iO;

    if (itype == typeO) {
      iO  = i;
      iH1 = hneigh[i*4  ];
      iH2 = hneigh[i*4+1];
      x1 = m[iO];
    }
    if (itype == typeH) {
      iO  = hneigh[i *4  ];
      iH1 = hneigh[iO*4  ];
      iH2 = hneigh[iO*4+1];
      if (iO >= inum) {
        non_local_oxy = 1;
      }
    }

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj,factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype = jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype = itype*lj_types+jtype;
      if (rsq < lj1[mtype].z) { // cut_ljsq
        numtyp r2inv = ucl_recip(rsq);
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp forcelj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        forcelj *= r2inv*factor_lj;

        f.x += delx*forcelj;
        f.y += dely*forcelj;
        f.z += delz*forcelj;

        if (EVFLAG && eflag) {
          numtyp e = r6inv * (lj3[mtype].x*r6inv-lj3[mtype].y);
          energy += factor_lj * (e - lj3[mtype].z);
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*forcelj;
          virial[1] += dely*dely*forcelj;
          virial[2] += delz*delz*forcelj;
          virial[3] += delx*dely*forcelj;
          virial[4] += delx*delz*forcelj;
          virial[5] += dely*delz*forcelj;
        }
      } // if LJ

      if (rsq < cut_coulsqplus) { //cut_coulsqplus
        int jH1, jH2, jO;
        numtyp qj; fetch(qj,j,q_tex);
        numtyp4 x2 = jx;
        if(itype == typeO || jtype == typeO) {
          if (jtype == typeO) {
            jO = j;
            jH1 = hneigh[j*4  ];
            jH2 = hneigh[j*4+1];
            x2 = m[j];
          }
          delx = x1.x-x2.x;
          dely = x1.y-x2.y;
          delz = x1.z-x2.z;
          rsq = delx*delx+dely*dely+delz*delz;
        }
        if (rsq < cut_coulsq) {
          numtyp r2inv = ucl_recip(rsq);
          numtyp r = ucl_rsqrt(r2inv);
          numtyp grij = g_ewald * r;
          numtyp expm2 = ucl_exp(-grij*grij);
          numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
          numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          numtyp prefactor = qj;
          prefactor *= qqrd2e*qtmp/r;
          numtyp force_coul = r2inv*prefactor * (_erfc + EWALD_F*grij*expm2 - factor_coul);

          if (itype == typeO) {
            fO.x += delx * force_coul;
            fO.y += dely * force_coul;
            fO.z += delz * force_coul;
            fO.w += 0;
          } else {
            f.x += delx * force_coul;
            f.y += dely * force_coul;
            f.z += delz * force_coul;
            f.w += 0;
          }
          if (EVFLAG && eflag) {
            e_coul += prefactor*(_erfc-factor_coul);
          }
          if (EVFLAG && vflag) {
            acctyp4 fd;
            fd.x = delx*force_coul;
            fd.y = dely*force_coul;
            fd.z = delz*force_coul;
            if (itype == typeO) {
              numtyp cO = 1 - alpha, cH = 0.5*alpha;
              numtyp4 vdi, vdj;
              numtyp4 xH1; fetch4(xH1,iH1,pos_tex);
              numtyp4 xH2; fetch4(xH2,iH2,pos_tex);
              numtyp4 xO; fetch4(xO,iO,pos_tex);
              vdi.x = xO.x*cO + xH1.x*cH + xH2.x*cH;
              vdi.y = xO.y*cO + xH1.y*cH + xH2.y*cH;
              vdi.z = xO.z*cO + xH1.z*cH + xH2.z*cH;
              //vdi.w = vdi.w;
              if (jtype == typeO) {
                numtyp4 xjH1; fetch4(xjH1,jH1,pos_tex);
                numtyp4 xjH2; fetch4(xjH2,jH2,pos_tex);
                numtyp4 xjO; fetch4(xjO,jO,pos_tex);
                vdj.x = xjO.x*cO + xjH1.x*cH + xjH2.x*cH;
                vdj.y = xjO.y*cO + xjH1.y*cH + xjH2.y*cH;
                vdj.z = xjO.z*cO + xjH1.z*cH + xjH2.z*cH;
                //vdj.w = vdj.w;
              } else vdj = jx;
              vO[0] += 0.5*(vdi.x - vdj.x)*fd.x;
              vO[1] += 0.5*(vdi.y - vdj.y)*fd.y;
              vO[2] += 0.5*(vdi.z - vdj.z)*fd.z;
              vO[3] += 0.5*(vdi.x - vdj.x)*fd.y;
              vO[4] += 0.5*(vdi.x - vdj.x)*fd.z;
              vO[5] += 0.5*(vdi.y - vdj.y)*fd.z;
            } else {
              if (jtype == typeO) {
                numtyp cO = 1 - alpha, cH = 0.5*alpha;
                numtyp4 vdj;
                numtyp4 xjH1; fetch4(xjH1,jH1,pos_tex);
                numtyp4 xjH2; fetch4(xjH2,jH2,pos_tex);
                numtyp4 xjO; fetch4(xjO,jO,pos_tex);
                vdj.x = xjO.x*cO + xjH1.x*cH + xjH2.x*cH;
                vdj.y = xjO.y*cO + xjH1.y*cH + xjH2.y*cH;
                vdj.z = xjO.z*cO + xjH1.z*cH + xjH2.z*cH;
                //vdj.w = vdj.w;
                virial[0] += (ix.x - vdj.x)*fd.x;
                virial[1] += (ix.y - vdj.y)*fd.y;
                virial[2] += (ix.z - vdj.z)*fd.z;
                virial[3] += (ix.x - vdj.x)*fd.y;
                virial[4] += (ix.x - vdj.x)*fd.z;
                virial[5] += (ix.y - vdj.y)*fd.z;
              } else {
                virial[0] += delx*fd.x;
                virial[1] += dely*fd.y;
                virial[2] += delz*fd.z;
                virial[3] += delx*fd.y;
                virial[4] += delx*fd.z;
                virial[5] += dely*fd.z;
              }
            }
          }
        }
        if (non_local_oxy == 1) {
          if (iO == j) {
            x2 = ix;
            qj = qtmp;
          }
          numtyp4 x1m = m[iO];
          delx = x1m.x-x2.x;
          dely = x1m.y-x2.y;
          delz = x1m.z-x2.z;
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cut_coulsq) {
            numtyp r2inv = ucl_recip(rsq);
            numtyp r = ucl_rsqrt(r2inv);
            numtyp grij = g_ewald * r;
            numtyp expm2 = ucl_exp(-grij*grij);
            numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
            numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

            numtyp prefactor = qj;
            prefactor *= qqrd2e*x1m.w/r;
            numtyp force_coul = r2inv*prefactor * (_erfc + EWALD_F*grij*expm2 - factor_coul);

            numtyp cO = 1 - alpha, cH = 0.5*alpha;
            numtyp4 fd;
            fd.x = delx * force_coul * cH;
            fd.y = dely * force_coul * cH;
            fd.z = delz * force_coul * cH;

            f.x += fd.x;
            f.y += fd.y;
            f.z += fd.z;

            if (EVFLAG && eflag) {
              e_coul += prefactor*(_erfc-factor_coul) * (acctyp)0.5 * alpha;
            }
            if (EVFLAG && vflag) {
              numtyp4 xH1; fetch4(xH1,iH1,pos_tex);
              numtyp4 xH2; fetch4(xH2,iH2,pos_tex);
              numtyp4 xO;  fetch4(xO,iO,pos_tex);

              virial[0] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.x;
              virial[1] += ((xO.y*cO + xH1.y*cH + xH2.y*cH) - x2.y) * fd.y;
              virial[2] += ((xO.z*cO + xH1.z*cH + xH2.z*cH) - x2.z) * fd.z;
              virial[3] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.y;
              virial[4] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.z;
              virial[5] += ((xO.y*cO + xH1.y*cH + xH2.y*cH) - x2.y) * fd.z;
            }
          }
        }
      } // if cut_coulsqplus
    } // for nbor
  } // if ii
  if (t_per_atom>1) {
#if (SHUFFLE_AVAIL == 0)
    red_acc[0][tid]=fO.x;
    red_acc[1][tid]=fO.y;
    red_acc[2][tid]=fO.z;
    red_acc[3][tid]=fO.w;
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
      simdsync();
      if (offset < s) {
        for (int r=0; r<4; r++)
          red_acc[r][tid] += red_acc[r][tid+s];
      }
    }
    fO.x=red_acc[0][tid];
    fO.y=red_acc[1][tid];
    fO.z=red_acc[2][tid];
    fO.w=red_acc[3][tid];
    if (EVFLAG && vflag) {
      simdsync();
      for (int r=0; r<6; r++) red_acc[r][tid]=vO[r];
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
        simdsync();
        if (offset < s) {
          for (int r=0; r<6; r++)
            red_acc[r][tid] += red_acc[r][tid+s];
        }
      }
      for (int r=0; r<6; r++) vO[r]=red_acc[r][tid];
    }
#else
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
      fO.x += shfl_down(fO.x, s, t_per_atom);
      fO.y += shfl_down(fO.y, s, t_per_atom);
      fO.z += shfl_down(fO.z, s, t_per_atom);
      fO.w += shfl_down(fO.w, s, t_per_atom);
    }
    if (EVFLAG && vflag) {
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
        for (int r=0; r<6; r++)
          vO[r] += shfl_down(vO[r], s, t_per_atom);
      }
    }
#endif
  }
  if(offset == 0 && ii<inum) {
    ansO[i] = fO;
    if (EVFLAG && vflag) {
      ansO[inum   + i].x = vO[0];
      ansO[inum   + i].y = vO[1];
      ansO[inum   + i].z = vO[2];
      ansO[inum*2 + i].x = vO[3];
      ansO[inum*2 + i].y = vO[4];
      ansO[inum*2 + i].z = vO[5];
    }
  }
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

__kernel void k_lj_tip4p_long_fast(const __global numtyp4 *restrict x_,
    const __global numtyp4 *restrict lj1_in,
    const __global numtyp4 *restrict lj3_in,
    const int lj_types,
    const __global numtyp *restrict sp_lj_in,
    const __global int * dev_nbor,
    const __global int * dev_packed,
    __global acctyp4 *restrict ans,
    __global acctyp *restrict engv,
    const int eflag, const int vflag, const int inum,
    const int nbor_pitch, const int t_per_atom,
    __global int *restrict hneigh,
    __global numtyp4 *restrict m,
    const int typeO, const int typeH,
    const numtyp alpha,
    const __global numtyp *restrict q_,
    const __global numtyp *restrict cutsq,
    const numtyp qqrd2e, const numtyp g_ewald,
    const numtyp cut_coulsq, const numtyp cut_coulsqplus,
    __global acctyp4 *restrict ansO) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (EVFLAG && eflag)
      lj3[tid]=lj3_in[tid];
  }
  acctyp4 f, fO;
  f.x=(acctyp)0;  f.y=(acctyp)0;  f.z=(acctyp)0;
  fO.x=(acctyp)0; fO.y=(acctyp)0; fO.z=(acctyp)0;
  acctyp energy, e_coul, virial[6], vO[6];
  if (EVFLAG) {
    energy = (acctyp)0;
    e_coul = (acctyp)0;
    for (int i=0; i<6; i++) {
      virial[i]=(acctyp)0;
      vO[i]=(acctyp)0;
    }
  }

  __syncthreads();
   if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
        n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype = ix.w;
    numtyp4 x1 = ix;

    int non_local_oxy = 0;
    int iH1, iH2, iO;

    if (itype == typeO) {
      iO  = i;
      iH1 = hneigh[i*4  ];
      iH2 = hneigh[i*4+1];
      x1 = m[iO];
    }
    if (itype == typeH) {
      iO  = hneigh[i *4  ];
      iH1 = hneigh[iO*4  ];
      iH2 = hneigh[iO*4+1];
      if (iO >= inum) {
        non_local_oxy = 1;
      }
    }

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj,factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype = jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype = itype*lj_types+jtype;
      if (rsq < lj1[mtype].z) { // cut_ljsq
        numtyp r2inv = ucl_recip(rsq);
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp forcelj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        forcelj *= r2inv*factor_lj;

        f.x += delx*forcelj;
        f.y += dely*forcelj;
        f.z += delz*forcelj;

        if (EVFLAG && eflag) {
          numtyp e = r6inv * (lj3[mtype].x*r6inv-lj3[mtype].y);
          energy += factor_lj * (e - lj3[mtype].z);
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*forcelj;
          virial[1] += dely*dely*forcelj;
          virial[2] += delz*delz*forcelj;
          virial[3] += delx*dely*forcelj;
          virial[4] += delx*delz*forcelj;
          virial[5] += dely*delz*forcelj;
        }
      } // if LJ

      if (rsq < cut_coulsqplus) { //cut_coulsqplus
        int jH1, jH2, jO;
        numtyp qj; fetch(qj,j,q_tex);
        numtyp4 x2 = jx;
        if(itype == typeO || jtype == typeO) {
          if (jtype == typeO) {
            jO = j;
            jH1 = hneigh[j*4  ];
            jH2 = hneigh[j*4+1];
            x2 = m[j];
          }
          delx = x1.x-x2.x;
          dely = x1.y-x2.y;
          delz = x1.z-x2.z;
          rsq = delx*delx+dely*dely+delz*delz;
        }
        if (rsq < cut_coulsq) {
          numtyp r2inv = ucl_recip(rsq);
          numtyp r = ucl_rsqrt(r2inv);
          numtyp grij = g_ewald * r;
          numtyp expm2 = ucl_exp(-grij*grij);
          numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
          numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          numtyp prefactor = qj;
          prefactor *= qqrd2e*qtmp/r;
          numtyp force_coul = r2inv*prefactor * (_erfc + EWALD_F*grij*expm2 - factor_coul);

          if (itype == typeO) {
            fO.x += delx * force_coul;
            fO.y += dely * force_coul;
            fO.z += delz * force_coul;
            fO.w += 0;
          } else {
            f.x += delx * force_coul;
            f.y += dely * force_coul;
            f.z += delz * force_coul;
            f.w += 0;
          }
          if (EVFLAG && eflag) {
            e_coul += prefactor*(_erfc-factor_coul);
          }
          if (EVFLAG && vflag) {
            acctyp4 fd;
            fd.x = delx*force_coul;
            fd.y = dely*force_coul;
            fd.z = delz*force_coul;
            if (itype == typeO) {
              numtyp cO = 1 - alpha, cH = 0.5*alpha;
              numtyp4 vdi, vdj;
              numtyp4 xH1; fetch4(xH1,iH1,pos_tex);
              numtyp4 xH2; fetch4(xH2,iH2,pos_tex);
              numtyp4 xO; fetch4(xO,iO,pos_tex);
              vdi.x = xO.x*cO + xH1.x*cH + xH2.x*cH;
              vdi.y = xO.y*cO + xH1.y*cH + xH2.y*cH;
              vdi.z = xO.z*cO + xH1.z*cH + xH2.z*cH;
              //vdi.w = vdi.w;
              if (jtype == typeO) {
                numtyp4 xjH1; fetch4(xjH1,jH1,pos_tex);
                numtyp4 xjH2; fetch4(xjH2,jH2,pos_tex);
                numtyp4 xjO; fetch4(xjO,jO,pos_tex);
                vdj.x = xjO.x*cO + xjH1.x*cH + xjH2.x*cH;
                vdj.y = xjO.y*cO + xjH1.y*cH + xjH2.y*cH;
                vdj.z = xjO.z*cO + xjH1.z*cH + xjH2.z*cH;
                //vdj.w = vdj.w;
              } else vdj = jx;
              vO[0] += 0.5*(vdi.x - vdj.x)*fd.x;
              vO[1] += 0.5*(vdi.y - vdj.y)*fd.y;
              vO[2] += 0.5*(vdi.z - vdj.z)*fd.z;
              vO[3] += 0.5*(vdi.x - vdj.x)*fd.y;
              vO[4] += 0.5*(vdi.x - vdj.x)*fd.z;
              vO[5] += 0.5*(vdi.y - vdj.y)*fd.z;
            } else {
              if (jtype == typeO) {
                numtyp cO = 1 - alpha, cH = 0.5*alpha;
                numtyp4 vdj;
                numtyp4 xjH1; fetch4(xjH1,jH1,pos_tex);
                numtyp4 xjH2; fetch4(xjH2,jH2,pos_tex);
                numtyp4 xjO; fetch4(xjO,jO,pos_tex);
                vdj.x = xjO.x*cO + xjH1.x*cH + xjH2.x*cH;
                vdj.y = xjO.y*cO + xjH1.y*cH + xjH2.y*cH;
                vdj.z = xjO.z*cO + xjH1.z*cH + xjH2.z*cH;
                //vdj.w = vdj.w;
                virial[0] += (ix.x - vdj.x)*fd.x;
                virial[1] += (ix.y - vdj.y)*fd.y;
                virial[2] += (ix.z - vdj.z)*fd.z;
                virial[3] += (ix.x - vdj.x)*fd.y;
                virial[4] += (ix.x - vdj.x)*fd.z;
                virial[5] += (ix.y - vdj.y)*fd.z;
              } else {
                virial[0] += delx*fd.x;
                virial[1] += dely*fd.y;
                virial[2] += delz*fd.z;
                virial[3] += delx*fd.y;
                virial[4] += delx*fd.z;
                virial[5] += dely*fd.z;
              }
            }
          }
        }
        if (non_local_oxy == 1) {
          if (iO == j) {
            x2 = ix;
            qj = qtmp;
          }
          numtyp4 x1m = m[iO];
          delx = x1m.x-x2.x;
          dely = x1m.y-x2.y;
          delz = x1m.z-x2.z;
          rsq = delx*delx+dely*dely+delz*delz;
          if (rsq < cut_coulsq) {
            numtyp r2inv = ucl_recip(rsq);
            numtyp r = ucl_rsqrt(r2inv);
            numtyp grij = g_ewald * r;
            numtyp expm2 = ucl_exp(-grij*grij);
            numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
            numtyp _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

            numtyp prefactor = qj;
            prefactor *= qqrd2e*x1m.w/r;
            numtyp force_coul = r2inv*prefactor * (_erfc + EWALD_F*grij*expm2 - factor_coul);

            numtyp cO = 1 - alpha, cH = 0.5*alpha;
            numtyp4 fd;
            fd.x = delx * force_coul * cH;
            fd.y = dely * force_coul * cH;
            fd.z = delz * force_coul * cH;

            f.x += fd.x;
            f.y += fd.y;
            f.z += fd.z;

            if (EVFLAG && eflag) {
              e_coul += prefactor*(_erfc-factor_coul) * (acctyp)0.5 * alpha;
            }
            if (EVFLAG && vflag) {
              numtyp4 xH1; fetch4(xH1,iH1,pos_tex);
              numtyp4 xH2; fetch4(xH2,iH2,pos_tex);
              numtyp4 xO;  fetch4(xO,iO,pos_tex);

              virial[0] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.x;
              virial[1] += ((xO.y*cO + xH1.y*cH + xH2.y*cH) - x2.y) * fd.y;
              virial[2] += ((xO.z*cO + xH1.z*cH + xH2.z*cH) - x2.z) * fd.z;
              virial[3] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.y;
              virial[4] += ((xO.x*cO + xH1.x*cH + xH2.x*cH) - x2.x) * fd.z;
              virial[5] += ((xO.y*cO + xH1.y*cH + xH2.y*cH) - x2.y) * fd.z;
            }
          }
        }
      } // if cut_coulsqplus
    } // for nbor
    if (t_per_atom>1) {
#if (SHUFFLE_AVAIL == 0)
      red_acc[0][tid]=fO.x;
      red_acc[1][tid]=fO.y;
      red_acc[2][tid]=fO.z;
      red_acc[3][tid]=fO.w;
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
        simdsync();
        if (offset < s) {
          for (int r=0; r<4; r++)
            red_acc[r][tid] += red_acc[r][tid+s];
        }
      }
      fO.x=red_acc[0][tid];
      fO.y=red_acc[1][tid];
      fO.z=red_acc[2][tid];
      fO.w=red_acc[3][tid];
      if (EVFLAG && vflag) {
        for (int r=0; r<6; r++) red_acc[r][tid]=vO[r];
        for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
          simdsync();
          if (offset < s) {
            for (int r=0; r<6; r++)
              red_acc[r][tid] += red_acc[r][tid+s];
          }
        }
        for (int r=0; r<6; r++) vO[r]=red_acc[r][tid];
      }
#else
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
        fO.x += shfl_down(fO.x, s, t_per_atom);
        fO.y += shfl_down(fO.y, s, t_per_atom);
        fO.z += shfl_down(fO.z, s, t_per_atom);
        fO.w += shfl_down(fO.w, s, t_per_atom);
      }
      if (EVFLAG && vflag) {
        for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
          for (int r=0; r<6; r++)
            vO[r] += shfl_down(vO[r], s, t_per_atom);
        }
      }
#endif
    }
    if(offset == 0) {
      ansO[i] = fO;
      if (EVFLAG && vflag) {
        ansO[inum   + i].x = vO[0];
        ansO[inum   + i].y = vO[1];
        ansO[inum   + i].z = vO[2];
        ansO[inum*2 + i].x = vO[3];
        ansO[inum*2 + i].y = vO[4];
        ansO[inum*2 + i].z = vO[5];
      }
    }
  } // if ii
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}
