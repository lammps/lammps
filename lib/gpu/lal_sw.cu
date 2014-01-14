// **************************************************************************
//                                   sw.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for acceleration of the sw pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Tue March 26, 2013
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"

#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float4> sw1_tex;
texture<float4> sw2_tex;
texture<float4> sw3_tex;
#else
texture<int4,1> pos_tex;
texture<int4> sw1_tex;
texture<int4> sw2_tex;
texture<int4> sw3_tex;
#endif

#else
#define pos_tex x_
#define sw1_tex sw1
#define sw2_tex sw2
#define sw3_tex sw3
#endif

#define THIRD (numtyp)0.66666667

#if (ARCH < 300)

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom, offset, \
                      eflag, vflag, ans, engv)                              \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_ELLIPSE];                               \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<4; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
    if (vflag>0) {                                                          \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<6; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]+=virial[i]*(acctyp)0.5;                                    \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#else

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom, offset, \
                      eflag, vflag, ans, engv)                              \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        energy += shfl_xor(energy, s, t_per_atom);                          \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]+=virial[i]*(acctyp)0.5;                                    \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#endif


__kernel void k_sw(const __global numtyp4 *restrict x_,
                   const __global numtyp4 *restrict sw1,
                   const __global numtyp4 *restrict sw2,
                   const __global numtyp4 *restrict sw3,
                   const __global int *restrict map,
                   const __global int *restrict elem2param,
                   const int nelements,
                   const __global int * dev_nbor, 
                   const __global int * dev_packed, 
                   __global acctyp4 *restrict ans, 
                   __global acctyp *restrict engv, 
                   const int eflag, const int vflag, const int inum, 
                   const int nbor_pitch, const int t_per_atom) {
  __local int n_stride;
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];
    
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      if (rsq<sw3[ijparam].y) { // sw_cutsq = sw3[ijparam].y
        numtyp4 sw1_ijparam; fetch4(sw1_ijparam,ijparam,sw1_tex);
        numtyp sw_epsilon=sw1_ijparam.x;
        numtyp sw_sigma=sw1_ijparam.y;
        numtyp4 sw2_ijparam; fetch4(sw2_ijparam,ijparam,sw2_tex);
        numtyp sw_biga=sw2_ijparam.x;
        numtyp sw_bigb=sw2_ijparam.y;
        numtyp sw_powerp=sw2_ijparam.z;
        numtyp sw_powerq=sw2_ijparam.w;
        numtyp4 sw3_ijparam; fetch4(sw3_ijparam,ijparam,sw3_tex);
        numtyp sw_cut=sw3_ijparam.x;
        numtyp sw_cutsq=sw3_ijparam.y;
        numtyp pre_sw_c1=sw_biga*sw_epsilon*sw_powerp*sw_bigb*
            pow(sw_sigma,sw_powerp);
        numtyp pre_sw_c2=sw_biga*sw_epsilon*sw_powerq*
            pow(sw_sigma,sw_powerq);
        numtyp pre_sw_c3=sw_biga*sw_epsilon*sw_bigb*
            pow(sw_sigma,sw_powerp+(numtyp)1.0);
        numtyp pre_sw_c4=sw_biga*sw_epsilon*
            pow(sw_sigma,sw_powerq+(numtyp)1.0);
        numtyp pre_sw_c5=sw_biga*sw_epsilon*sw_bigb*
            pow(sw_sigma,sw_powerp);
        numtyp pre_sw_c6=sw_biga*sw_epsilon*
            pow(sw_sigma,sw_powerq);

        numtyp r=ucl_sqrt(rsq);
        numtyp rp=ucl_powr(r,-sw_powerp);
        numtyp rq=ucl_powr(r,-sw_powerq);
        numtyp rainv=ucl_recip(r-sw_cut);
        numtyp expsrainv=ucl_exp(sw_sigma*rainv);
        rainv*=rainv*r;
        numtyp force = (pre_sw_c1*rp-pre_sw_c2*rq +
                       (pre_sw_c3*rp-pre_sw_c4*rq) * rainv)*
                       expsrainv*ucl_recip(rsq);
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) 
          energy+=(pre_sw_c5*rp - pre_sw_c6*rq) * expsrainv; 

        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }
    } // for nbor

    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii

}

#define threebody(delr1x, delr1y, delr1z, eflag, energy)                     \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - sw_cut_ij);                                 \
  numtyp gsrainv1 = sw_sigma_gamma_ij * rainv1;                              \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rinvsq2 = ucl_recip(rsq2);                                          \
  numtyp rainv2 = ucl_recip(r2 - sw_cut_ik);                                 \
  numtyp gsrainv2 = sw_sigma_gamma_ik * rainv2;                              \
  numtyp gsrainvsq2 = gsrainv2*rainv2/r2;                                    \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - sw_costheta_ijk;                                       \
  numtyp delcssq = delcs*delcs;                                              \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = sw_lambda_epsilon_ijk * facexp*delcssq;                    \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp frad2 = facrad*gsrainvsq2;                                          \
  numtyp facang = sw_lambda_epsilon2_ijk * facexp*delcs;                     \
  numtyp facang12 = rinv12*facang;                                           \
  numtyp csfacang = cs*facang;                                               \
  numtyp csfac1 = rinvsq1*csfacang;                                          \
                                                                             \
  fjx = delr1x*(frad1+csfac1)-delr2x*facang12;                               \
  fjy = delr1y*(frad1+csfac1)-delr2y*facang12;                               \
  fjz = delr1z*(frad1+csfac1)-delr2z*facang12;                               \
                                                                             \
  numtyp csfac2 = rinvsq2*csfacang;                                          \
                                                                             \
  fkx = delr2x*(frad2+csfac2)-delr1x*facang12;                               \
  fky = delr2y*(frad2+csfac2)-delr1y*facang12;                               \
  fkz = delr2z*(frad2+csfac2)-delr1z*facang12;                               \
                                                                             \
  if (eflag>0)                                                               \
    energy+=facrad;                                                          \
  if (vflag>0) {                                                             \
    virial[0] += delr1x*fjx + delr2x*fkx;                                    \
    virial[1] += delr1y*fjy + delr2y*fky;                                    \
    virial[2] += delr1z*fjz + delr2z*fkz;                                    \
    virial[3] += delr1x*fjy + delr2x*fky;                                    \
    virial[4] += delr1x*fjz + delr2x*fkz;                                    \
    virial[5] += delr1y*fjz + delr2y*fkz;                                    \
  }                                                                          \
}

#define threebody_half(delr1x, delr1y, delr1z)                               \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - sw_cut_ij);                                 \
  numtyp gsrainv1 = sw_sigma_gamma_ij * rainv1;                              \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rainv2 = ucl_recip(r2 - sw_cut_ik);                                 \
  numtyp gsrainv2 = sw_sigma_gamma_ik * rainv2;                              \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - sw_costheta_ijk;                                       \
  numtyp delcssq = delcs*delcs;                                              \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = sw_lambda_epsilon_ijk * facexp*delcssq;                    \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp facang = sw_lambda_epsilon2_ijk * facexp*delcs;                     \
  numtyp facang12 = rinv12*facang;                                           \
  numtyp csfacang = cs*facang;                                               \
  numtyp csfac1 = rinvsq1*csfacang;                                          \
                                                                             \
  fjx = delr1x*(frad1+csfac1)-delr2x*facang12;                               \
  fjy = delr1y*(frad1+csfac1)-delr2y*facang12;                               \
  fjz = delr1z*(frad1+csfac1)-delr2z*facang12;                               \
}

__kernel void k_sw_three_center(const __global numtyp4 *restrict x_, 
                                const __global numtyp4 *restrict sw1,
                                const __global numtyp4 *restrict sw2,
                                const __global numtyp4 *restrict sw3,
                                const __global int *restrict map,
                                const __global int *restrict elem2param,
                                const int nelements,
                                const __global int * dev_nbor, 
                                const __global int * dev_packed, 
                                __global acctyp4 *restrict ans, 
                                __global acctyp *restrict engv, 
                                const int eflag, const int vflag, 
                                const int inum,  const int nbor_pitch, 
                                const int t_per_atom, const int evatom) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp sw_epsilon, sw_sigma, sw_lambda, sw_gamma;
  numtyp sw_sigma_gamma_ij, sw_cut_ij, sw_sigma_gamma_ik, sw_cut_ik;
  numtyp sw_costheta_ijk, sw_lambda_epsilon_ijk, sw_lambda_epsilon2_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();
  
  if (ii<inum) {
    const __global int *nbor_j, *list_end;
    int i, numj;

    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,list_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w; 
    itype=map[itype];

    for ( ; nbor_j<list_end; nbor_j+=n_stride) {
  
      int j=*nbor_j;
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = jx.x-ix.x;
      numtyp delr1y = jx.y-ix.y;
      numtyp delr1z = jx.z-ix.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];
      numtyp4 sw3_ijparam; fetch4(sw3_ijparam,ijparam,sw3_tex);

      if (rsq1 > sw3_ijparam.y) continue;

      numtyp4 sw1_ijparam; fetch4(sw1_ijparam,ijparam,sw1_tex);
      sw_sigma=sw1_ijparam.y;
      sw_gamma=sw1_ijparam.w;
      sw_sigma_gamma_ij=sw1_ijparam.y*sw1_ijparam.w; //sw_sigma*sw_gamma;
      sw_cut_ij=sw3_ijparam.x;

      const __global int *nbor_k=nbor_j-offset_j+offset_k;
      if (nbor_k<=nbor_j)
        nbor_k+=n_stride;

      for ( ; nbor_k<list_end; nbor_k+=n_stride) {
        int k=*nbor_k;
        k &= NEIGHMASK;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[itype*nelements*nelements+ktype*nelements+ktype];
        numtyp4 sw3_ikparam; fetch4(sw3_ikparam,ikparam,sw3_tex);

        numtyp delr2x = kx.x-ix.x;
        numtyp delr2y = kx.y-ix.y;
        numtyp delr2z = kx.z-ix.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;
        if (rsq2 < sw3_ikparam.y) {   // sw_cutsq=sw3[ikparam].y;
          numtyp4 sw1_ikparam; fetch4(sw1_ikparam,ikparam,sw1_tex);
          sw_sigma=sw1_ikparam.y;
          sw_gamma=sw1_ikparam.w;
          sw_sigma_gamma_ik=sw1_ikparam.y*sw1_ikparam.w; //sw_sigma*sw_gamma;
          sw_cut_ik=sw3_ikparam.x;

          int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];
          numtyp4 sw1_ijkparam; fetch4(sw1_ijkparam,ijkparam,sw1_tex);
          sw_epsilon=sw1_ijkparam.x;
          sw_lambda=sw1_ijkparam.z;
          sw_lambda_epsilon_ijk=sw1_ijkparam.x*sw1_ijkparam.z; //sw_lambda*sw_epsilon;
          sw_lambda_epsilon2_ijk=(numtyp)2.0*sw_lambda_epsilon_ijk;
          numtyp4 sw3_ijkparam; fetch4(sw3_ijkparam,ijkparam,sw3_tex);
          sw_costheta_ijk=sw3_ijkparam.z;

          numtyp fjx, fjy, fjz, fkx, fky, fkz;
          threebody(delr1x,delr1y,delr1z,eflag,energy);

          f.x -= fjx + fkx;
          f.y -= fjy + fky;
          f.z -= fjz + fkz;
        }  
      }
    } // for nbor

    numtyp pre;
    if (evatom==1)
      pre=THIRD;
    else
      pre=(numtyp)2.0;
    energy*=pre;
    for (int i=0; i<6; i++)
      virial[i]*=pre;

    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                    eflag,vflag,ans,engv);

  } // if ii
}

__kernel void k_sw_three_end(const __global numtyp4 *restrict x_, 
                             const __global numtyp4 *restrict sw1,
                             const __global numtyp4 *restrict sw2,
                             const __global numtyp4 *restrict sw3,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements,
                             const __global int * dev_nbor, 
                             const __global int * dev_packed, 
                             __global acctyp4 *restrict ans, 
                             __global acctyp *restrict engv, 
                             const int eflag, const int vflag, 
                             const int inum,  const int nbor_pitch, 
                             const int t_per_atom) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp sw_epsilon, sw_sigma, sw_lambda, sw_gamma;
  numtyp sw_sigma_gamma_ij, sw_cut_ij, sw_sigma_gamma_ik, sw_cut_ik;
  numtyp sw_costheta_ijk, sw_lambda_epsilon_ijk, sw_lambda_epsilon2_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();
  
  if (ii<inum) {
    const __global int *nbor_j, *list_end, *k_end;
    int i, numj;

    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,list_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    for ( ; nbor_j<list_end; nbor_j+=n_stride) {
      int j=*nbor_j;
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];
      numtyp4 sw3_ijparam; fetch4(sw3_ijparam,ijparam,sw3_tex);

      if (rsq1 > sw3_ijparam.y) continue;

      numtyp4 sw1_ijparam; fetch4(sw1_ijparam,ijparam,sw1_tex);
      sw_sigma=sw1_ijparam.y;
      sw_gamma=sw1_ijparam.w;
      sw_sigma_gamma_ij=sw1_ijparam.y*sw1_ijparam.w; //sw_sigma*sw_gamma;
      sw_cut_ij=sw3_ijparam.x;

      const __global int *nbor_k=dev_nbor+j+nbor_pitch;
      int numk=*nbor_k; 
      if (dev_nbor==dev_packed) {
        nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
        k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
        nbor_k+=offset_k;
      } else {
        nbor_k+=nbor_pitch;
        nbor_k=dev_packed+*nbor_k;
        k_end=nbor_k+numk;
        nbor_k+=offset_k;
      }

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=*nbor_k;
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[itype*nelements*nelements+ktype*nelements+ktype]; 

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;
        numtyp4 sw3_ikparam; fetch4(sw3_ikparam,ikparam,sw3_tex);

        if (rsq2 < sw3_ikparam.y) {
          numtyp4 sw1_ikparam; fetch4(sw1_ikparam,ikparam,sw1_tex);
          sw_sigma=sw1_ikparam.y;
          sw_gamma=sw1_ikparam.w;
          sw_sigma_gamma_ik=sw1_ikparam.y*sw1_ikparam.w; //sw_sigma*sw_gamma;
          sw_cut_ik=sw3_ikparam.x;

          int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];
          numtyp4 sw1_ijkparam; fetch4(sw1_ijkparam,ijkparam,sw1_tex);
          sw_epsilon=sw1_ijkparam.x;
          sw_lambda=sw1_ijkparam.z;
          sw_lambda_epsilon_ijk=sw1_ijkparam.x*sw1_ijkparam.z; //sw_lambda*sw_epsilon;
          sw_lambda_epsilon2_ijk=(numtyp)2.0*sw_lambda_epsilon_ijk;
          numtyp4 sw3_ijkparam; fetch4(sw3_ijkparam,ijkparam,sw3_tex);
          sw_costheta_ijk=sw3_ijkparam.z;

          numtyp fjx, fjy, fjz;
          //if (evatom==0) {
            threebody_half(delr1x,delr1y,delr1z);
          //} else {
          //  numtyp fkx, fky, fkz;
          //  threebody(delr1x,delr1y,delr1z,eflag,energy);
          //}

          f.x += fjx;
          f.y += fjy;
          f.z += fjz;
        }
      }

    } // for nbor
    #ifdef THREE_CONCURRENT
    store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv);
    #else
    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                    eflag,vflag,ans,engv);
    #endif
  } // if ii
}

__kernel void k_sw_three_end_vatom(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict sw1,
                             const __global numtyp4 *restrict sw2,
                             const __global numtyp4 *restrict sw3,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements,
                             const __global int * dev_nbor, 
                             const __global int * dev_packed, 
                             __global acctyp4 *restrict ans, 
                             __global acctyp *restrict engv, 
                             const int eflag, const int vflag, 
                             const int inum,  const int nbor_pitch, 
                             const int t_per_atom) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp sw_epsilon, sw_sigma, sw_lambda, sw_gamma;
  numtyp sw_sigma_gamma_ij, sw_cut_ij, sw_sigma_gamma_ik, sw_cut_ik;
  numtyp sw_costheta_ijk, sw_lambda_epsilon_ijk, sw_lambda_epsilon2_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();
  
  if (ii<inum) {
    const __global int *nbor_j, *list_end, *k_end;
    int i, numj;

    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,list_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    for ( ; nbor_j<list_end; nbor_j+=n_stride) {
      int j=*nbor_j;
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];
      numtyp4 sw3_ijparam; fetch4(sw3_ijparam,ijparam,sw3_tex);

      if (rsq1 > sw3_ijparam.y) continue;

      numtyp4 sw1_ijparam; fetch4(sw1_ijparam,ijparam,sw1_tex);
      sw_sigma=sw1_ijparam.y;
      sw_gamma=sw1_ijparam.w;
      sw_sigma_gamma_ij=sw1_ijparam.y*sw1_ijparam.w; //sw_sigma*sw_gamma;
      sw_cut_ij=sw3_ijparam.x;
        
      const __global int *nbor_k=dev_nbor+j+nbor_pitch;
      int numk=*nbor_k;
      if (dev_nbor==dev_packed) {
        nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
        k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
        nbor_k+=offset_k;
      } else {
        nbor_k+=nbor_pitch;
        nbor_k=dev_packed+*nbor_k;
        k_end=nbor_k+numk;
        nbor_k+=offset_k;
      }

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=*nbor_k;
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[itype*nelements*nelements+ktype*nelements+ktype]; 
        numtyp4 sw3_ikparam; fetch4(sw3_ikparam,ikparam,sw3_tex);

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;

        if (rsq2 < sw3_ikparam.y) {
          numtyp4 sw1_ikparam; fetch4(sw1_ikparam,ikparam,sw1_tex);
          sw_sigma=sw1_ikparam.y;
          sw_gamma=sw1_ikparam.w;
          sw_sigma_gamma_ik=sw1_ikparam.y*sw1_ikparam.w; //sw_sigma*sw_gamma;
          sw_cut_ik=sw3_ikparam.x;

          int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];
          numtyp4 sw1_ijkparam; fetch4(sw1_ijkparam,ijkparam,sw1_tex);
          sw_epsilon=sw1_ijkparam.x;
          sw_lambda=sw1_ijkparam.z;
          sw_lambda_epsilon_ijk=sw1_ijkparam.x*sw1_ijkparam.z; //sw_lambda*sw_epsilon;
          sw_lambda_epsilon2_ijk=(numtyp)2.0*sw_lambda_epsilon_ijk;
          numtyp4 sw3_ijkparam; fetch4(sw3_ijkparam,ijkparam,sw3_tex);
          sw_costheta_ijk=sw3_ijkparam.z;

          numtyp fjx, fjy, fjz, fkx, fky, fkz;
          threebody(delr1x,delr1y,delr1z,eflag,energy);

          f.x += fjx;
          f.y += fjy;
          f.z += fjz;
        }
      }

    } // for nbor
    energy*=THIRD;
    for (int i=0; i<6; i++)
      virial[i]*=THIRD;
    #ifdef THREE_CONCURRENT
    store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv);
    #else
    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                    eflag,vflag,ans,engv);
    #endif
  } // if ii
}

