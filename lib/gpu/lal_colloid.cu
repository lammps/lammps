// **************************************************************************
//                                   colloid.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the colloid pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : nguyentd@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
#else
_texture_2d( pos_tex,int4);
#endif
#else
#define pos_tex x_
#endif

__kernel void k_colloid(const __global numtyp4 *restrict x_,
                        const __global numtyp4 *restrict lj1,
                        const __global numtyp4 *restrict lj3,
                        const int lj_types,
                        const __global numtyp *restrict sp_lj_in,
                        const __global numtyp4 *restrict colloid1,
                        const __global numtyp4 *restrict colloid2,
                        const __global int *form,
                        const __global int *dev_nbor,
                        const __global int *dev_packed,
                        __global acctyp4 *restrict ans,
                        __global acctyp *restrict engv,
                        const int eflag, const int vflag, const int inum,
                        const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[4];
  int n_stride;
  local_allocate_store_pair();

  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<lj1[mtype].z) {
        numtyp r,r2inv,r6inv;
        numtyp c1,c2,fR,evdwl;
        numtyp K[9],h[4],g[4];
        numtyp force = (numtyp)0;

        if (form[mtype]==0) { // SMALL_SMALL
          r2inv=ucl_recip(rsq);
          r6inv = r2inv*r2inv*r2inv;
          force = r2inv*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
          force*=factor_lj;
        } else if (form[mtype]==1) { // SMALL_LARGE
          c2 = colloid1[mtype].z;
          K[1] = c2*c2;
          K[2] = rsq;
          K[0] = K[1] - rsq;
          K[4] = rsq*rsq;
          K[3] = K[1] - K[2];
          K[3] *= K[3]*K[3];
          K[6] = K[3]*K[3];
          fR = colloid2[mtype].z*colloid1[mtype].x*c2*K[1]/K[3];
          force = (numtyp)4.0/(numtyp)15.0*fR *
             ((numtyp)2.0*(K[1]+K[2]) *
             (K[1]*((numtyp)5.0*K[1]+(numtyp)22.0*K[2])+(numtyp)5.0*K[4]) *
             colloid2[mtype].w/K[6]-(numtyp)5.0) / K[0];
          force*=factor_lj;
        } else if (form[mtype]==2) { // LARGE_LARGE
          r = ucl_sqrt(rsq);
          c1 = colloid1[mtype].y;
          c2 = colloid1[mtype].z;
          K[0] = c1*c2;
          K[1] = c1+c2;
          K[2] = c1-c2;
          K[3] = K[1]+r;
          K[4] = K[1]-r;
          K[5] = K[2]+r;
          K[6] = K[2]-r;
          K[7] = ucl_recip(K[3]*K[4]);
          K[8] = ucl_recip(K[5]*K[6]);
          g[0] = (numtyp)1.0/(K[3]*K[3]*K[3]*K[3]*K[3]*K[3]*K[3]);  // ucl_powr(K[3],(numtyp)-7.0);
          g[1] = (numtyp)1.0/(K[4]*K[4]*K[4]*K[4]*K[4]*K[4]*K[4]);  //-ucl_powr(-K[4],(numtyp)-7.0);
          g[2] = (numtyp)1.0/(K[5]*K[5]*K[5]*K[5]*K[5]*K[5]*K[5]);  // ucl_powr(K[5],(numtyp)-7.0);
          g[3] = (numtyp)1.0/(K[6]*K[6]*K[6]*K[6]*K[6]*K[6]*K[6]);  //-ucl_powr(-K[6],(numtyp)-7.0);
          h[0] = ((K[3]+(numtyp)5.0*K[1])*K[3]+(numtyp)30.0*K[0])*g[0];
          h[1] = ((K[4]+(numtyp)5.0*K[1])*K[4]+(numtyp)30.0*K[0])*g[1];
          h[2] = ((K[5]+(numtyp)5.0*K[2])*K[5]-(numtyp)30.0*K[0])*g[2];
          h[3] = ((K[6]+(numtyp)5.0*K[2])*K[6]-(numtyp)30.0*K[0])*g[3];
          g[0] *= (numtyp)42.0*K[0]/K[3]+(numtyp)6.0*K[1]+K[3];
          g[1] *= (numtyp)42.0*K[0]/K[4]+(numtyp)6.0*K[1]+K[4];
          g[2] *= (numtyp)-42.0*K[0]/K[5]+(numtyp)6.0*K[2]+K[5];
          g[3] *= (numtyp)-42.0*K[0]/K[6]+(numtyp)6.0*K[2]+K[6];

          fR = colloid1[mtype].x*colloid2[mtype].w/r/(numtyp)37800.0;
          evdwl = fR * (h[0]-h[1]-h[2]+h[3]);
          numtyp dUR = evdwl/r + (numtyp)5.0*fR*(g[0]+g[1]-g[2]-g[3]);
          numtyp dUA = -colloid1[mtype].x/(numtyp)3.0*r*
                       (((numtyp)2.0*K[0]*K[7]+(numtyp)1.0)*K[7] +
                       ((numtyp)2.0*K[0]*K[8]-(numtyp)1.0)*K[8]);
          force = factor_lj * (dUR+dUA)/r;
        }

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          numtyp e=(numtyp)0.0;
          if (form[mtype]==0) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          } else if (form[mtype]==1) {
            e=(numtyp)2.0/(numtyp)9.0*fR *
              ((numtyp)1.0-(K[1]*(K[1]*(K[1]/(numtyp)3.0+(numtyp)3.0*K[2]) +
              (numtyp)4.2*K[4])+K[2]*K[4]) * colloid2[mtype].w/K[6]);
          } else if (form[mtype]==2) {
            e=evdwl+colloid1[mtype].x/(numtyp)6.0 *
              ((numtyp)2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7]));
          }
          energy+=factor_lj*(e-lj3[mtype].z);
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}

__kernel void k_colloid_fast(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict lj1_in,
                             const __global numtyp4 *restrict lj3_in,
                             const __global numtyp *restrict sp_lj_in,
                             const __global numtyp4 *restrict colloid1_in,
                             const __global numtyp4 *restrict colloid2_in,
                             const __global int *form_in,
                             const __global int *dev_nbor,
                             const __global int *dev_packed,
                             __global acctyp4 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag, const int inum,
                             const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 colloid1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 colloid2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local int form[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  int n_stride;
  local_allocate_store_pair();

  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    colloid1[tid]=colloid1_in[tid];
    colloid2[tid]=colloid2_in[tid];
    form[tid]=form_in[tid];
    if (EVFLAG && eflag)
      lj3[tid]=lj3_in[tid];
  }

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<lj1[mtype].z) {
        numtyp r,r2inv,r6inv;
        numtyp c1,c2,fR,evdwl;
        numtyp K[9],h[4],g[4];
        numtyp force = (numtyp)0;

        if (form[mtype]==0) { // SMALL_SMALL
          r2inv=ucl_recip(rsq);
          r6inv = r2inv*r2inv*r2inv;
          force = r2inv*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
          force*=factor_lj;
        } else if (form[mtype]==1) { // SMALL_LARGE
          c2 = colloid1[mtype].z;
          K[1] = c2*c2;
          K[2] = rsq;
          K[0] = K[1] - rsq;
          K[4] = rsq*rsq;
          K[3] = K[1] - K[2];
          K[3] *= K[3]*K[3];
          K[6] = K[3]*K[3];
          fR = colloid2[mtype].z*colloid1[mtype].x*c2*K[1]/K[3];
          force = (numtyp)4.0/(numtyp)15.0*fR *
            ((numtyp)2.0*(K[1]+K[2]) *
            (K[1]*((numtyp)5.0*K[1]+(numtyp)22.0*K[2])+(numtyp)5.0*K[4]) *
            colloid2[mtype].w/K[6]-(numtyp)5.0) / K[0];
          force*=factor_lj;
        } else if (form[mtype]==2) { // LARGE_LARGE
          r = ucl_sqrt(rsq);
          c1 = colloid1[mtype].y;
          c2 = colloid1[mtype].z;
          K[0] = c1*c2;
          K[1] = c1+c2;
          K[2] = c1-c2;
          K[3] = K[1]+r;
          K[4] = K[1]-r;
          K[5] = K[2]+r;
          K[6] = K[2]-r;
          K[7] = ucl_recip(K[3]*K[4]);
          K[8] = ucl_recip(K[5]*K[6]);
          g[0] = (numtyp)1.0/(K[3]*K[3]*K[3]*K[3]*K[3]*K[3]*K[3]);  // ucl_powr(K[3],(numtyp)-7.0);
          g[1] = (numtyp)1.0/(K[4]*K[4]*K[4]*K[4]*K[4]*K[4]*K[4]);  //-ucl_powr(-K[4],(numtyp)-7.0);
          g[2] = (numtyp)1.0/(K[5]*K[5]*K[5]*K[5]*K[5]*K[5]*K[5]);  // ucl_powr(K[5],(numtyp)-7.0);
          g[3] = (numtyp)1.0/(K[6]*K[6]*K[6]*K[6]*K[6]*K[6]*K[6]);  //-ucl_powr(-K[6],(numtyp)-7.0);
          h[0] = ((K[3]+(numtyp)5.0*K[1])*K[3]+(numtyp)30.0*K[0])*g[0];
          h[1] = ((K[4]+(numtyp)5.0*K[1])*K[4]+(numtyp)30.0*K[0])*g[1];
          h[2] = ((K[5]+(numtyp)5.0*K[2])*K[5]-(numtyp)30.0*K[0])*g[2];
          h[3] = ((K[6]+(numtyp)5.0*K[2])*K[6]-(numtyp)30.0*K[0])*g[3];
          g[0] *= (numtyp)42.0*K[0]/K[3]+(numtyp)6.0*K[1]+K[3];
          g[1] *= (numtyp)42.0*K[0]/K[4]+(numtyp)6.0*K[1]+K[4];
          g[2] *= (numtyp)-42.0*K[0]/K[5]+(numtyp)6.0*K[2]+K[5];
          g[3] *= (numtyp)-42.0*K[0]/K[6]+(numtyp)6.0*K[2]+K[6];

          fR = colloid1[mtype].x*colloid2[mtype].w/r/(numtyp)37800.0;
          evdwl = fR * (h[0]-h[1]-h[2]+h[3]);
          numtyp dUR = evdwl/r + (numtyp)5.0*fR*(g[0]+g[1]-g[2]-g[3]);
          numtyp dUA = -colloid1[mtype].x/(numtyp)3.0*r*
            (((numtyp)2.0*K[0]*K[7]+(numtyp)1.0)*K[7] +
            ((numtyp)2.0*K[0]*K[8]-(numtyp)1.0)*K[8]);
          force = factor_lj * (dUR+dUA)/r;
        } else force = (numtyp)0.0;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          numtyp e=(numtyp)0.0;
          if (form[mtype]==0) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          } else if (form[mtype]==1) {
            e=(numtyp)2.0/(numtyp)9.0*fR *
              ((numtyp)1.0-(K[1]*(K[1]*(K[1]/(numtyp)3.0+
              (numtyp)3.0*K[2])+(numtyp)4.2*K[4])+K[2]*K[4])*
              colloid2[mtype].w/K[6]);
          } else if (form[mtype]==2) {
            e=evdwl+colloid1[mtype].x/(numtyp)6.0 *
              ((numtyp)2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7]));
          }
          energy+=factor_lj*(e-lj3[mtype].z);
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}

