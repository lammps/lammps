// **************************************************************************
//                                dipole_lj.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the dipole/cut pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float> q_tex;
texture<float4> mu_tex;
#else
texture<int4,1> pos_tex;
texture<int2> q_tex;
texture<int4,1> mu_tex;
#endif

#else
#define pos_tex x_
#define q_tex q_
#define mu_tex mu_
#endif

#if (ARCH < 300)

#define store_answers_tq(f, tor, energy, ecoul, virial, ii, inum, tid,      \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[8][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=tor.x;                                                  \
    red_acc[4][tid]=tor.y;                                                  \
    red_acc[5][tid]=tor.z;                                                  \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    tor.x=red_acc[3][tid];                                                  \
    tor.y=red_acc[4][tid];                                                  \
    tor.z=red_acc[5][tid];                                                  \
    if (eflag>0 || vflag>0) {                                               \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      red_acc[6][tid]=energy;                                               \
      red_acc[7][tid]=ecoul;                                                \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<8; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
      energy=red_acc[6][tid];                                               \
      ecoul=red_acc[7][tid];                                                \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    engv+=ii;                                                               \
    if (eflag>0) {                                                          \
      *engv=energy;                                                         \
      engv+=inum;                                                           \
      *engv=e_coul;                                                         \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *engv=virial[i];                                                    \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }

#else

#define store_answers_tq(f, tor, energy, e_coul, virial, ii, inum, tid,     \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        tor.x += shfl_xor(tor.x, s, t_per_atom);                            \
        tor.y += shfl_xor(tor.y, s, t_per_atom);                            \
        tor.z += shfl_xor(tor.z, s, t_per_atom);                            \
        energy += shfl_xor(energy, s, t_per_atom);                          \
        e_coul += shfl_xor(e_coul, s, t_per_atom);                          \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    engv+=ii;                                                               \
    if (eflag>0) {                                                          \
      *engv=energy;                                                         \
      engv+=inum;                                                           \
      *engv=e_coul;                                                         \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *engv=virial[i];                                                    \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }

#endif

__kernel void k_dipole_lj(const __global numtyp4 *restrict x_, 
                          const __global numtyp4 *restrict lj1,
                          const __global numtyp4 *restrict lj3, 
                          const int lj_types, 
                          const __global numtyp *restrict sp_lj_in, 
                          const __global int *dev_nbor, 
                          const __global int *dev_packed, 
                          __global acctyp4 *restrict ans,
                          __global acctyp *restrict engv, 
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch, 
                          const __global numtyp *restrict q_,
                          const __global numtyp4 *restrict mu_,
                          const __global numtyp *restrict cutsq,
                          const numtyp qqrd2e, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[8];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  sp_lj[4]=sp_lj_in[4];
  sp_lj[5]=sp_lj_in[5];
  sp_lj[6]=sp_lj_in[6];
  sp_lj[7]=sp_lj_in[7];

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp4 tor;
  tor.x=(acctyp)0;
  tor.y=(acctyp)0;
  tor.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    numtyp4 mui; fetch4(mui,i,mu_tex); //mu_[i];
    int itype=ix.w;

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      numtyp qj; fetch(qj,j,q_tex);
      numtyp4 muj; fetch4(muj,j,mu_tex); //mu_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, r6inv;
        numtyp rinv, r3inv, r5inv, r7inv;
        numtyp pre1, pre2, pre3, pre4;
        numtyp pdotp, pidotr, pjdotr;
        acctyp4 forcecoul, ticoul;
        acctyp4 force;

        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;
        } else force_lj = (numtyp)0.0;

        if (rsq < lj1[mtype].w) {
          rinv = ucl_rsqrt(rsq);

          // charge-charge
          if (qtmp != (numtyp)0.0 && qj != (numtyp)0.0) { 
            r3inv = r2inv*rinv;
            pre1 = qtmp*qj*r3inv;

            forcecoul.x += pre1*delx;
            forcecoul.y += pre1*dely;
            forcecoul.z += pre1*delz;
          }                    

          // dipole-dipole
          if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0) {
            r3inv = r2inv*rinv; 
            r5inv = r3inv*r2inv;
	          r7inv = r5inv*r2inv;
            pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

            pre1 = (numtyp)3.0*r5inv*pdotp - (numtyp)15.0*r7inv*pidotr*pjdotr;
            pre2 = (numtyp)3.0*r5inv*pjdotr;
            pre3 = (numtyp)3.0*r5inv*pidotr;
            pre4 = (numtyp)(-1.0)*r3inv;

            forcecoul.x += pre1*delx + pre2*mui.x + pre3*muj.x;
            forcecoul.y += pre1*dely + pre2*mui.y + pre3*muj.y;
            forcecoul.z += pre1*delz + pre2*mui.z + pre3*muj.z;
	    
            numtyp crossx = pre4 * (mui.y*muj.z - mui.z*muj.y);
            numtyp crossy = pre4 * (mui.z*muj.x - mui.x*muj.z);
            numtyp crossz = pre4 * (mui.x*muj.y - mui.y*muj.x);

            ticoul.x += crossx + pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += crossy + pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += crossz + pre2 * (mui.x*dely - mui.y*delx);
          }

          // dipole-charge
          if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) {
            r3inv = r2inv*rinv; 
            r5inv = r3inv*r2inv;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pre1 = (numtyp)3.0*qj*r5inv * pidotr;
            pre2 = qj*r3inv;
              
            forcecoul.x += pre2*mui.x - pre1*delx;
            forcecoul.y += pre2*mui.y - pre1*dely;
            forcecoul.z += pre2*mui.z - pre1*delz;
            ticoul.x += pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += pre2 * (mui.x*dely - mui.y*delx);
          }
          
          // charge-dipole
          if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;
            pre1 = (numtyp)3.0*qtmp*r5inv * pjdotr;
            pre2 = qtmp*r3inv;
            
            forcecoul.x += pre1*delx - pre2*muj.x;
            forcecoul.y += pre1*dely - pre2*muj.y;
            forcecoul.z += pre1*delz - pre2*muj.z;
          }
        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
          ticoul.x = ticoul.y = ticoul.z = (acctyp)0;
        }

        numtyp fq = factor_coul*qqrd2e;
        force.x = fq*forcecoul.x + delx*force_lj;
        force.y = fq*forcecoul.y + dely*force_lj;
        force.z = fq*forcecoul.z + delz*force_lj;
        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=fq*ticoul.x;
        tor.y+=fq*ticoul.y;
        tor.z+=fq*ticoul.z;

        if (eflag>0) {
          acctyp e = (acctyp)0.0;  
          if (rsq < lj1[mtype].w) {
            e = qtmp*qj*rinv;
            if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0)
              e += r3inv*pdotp - (numtyp)3.0*r5inv*pidotr*pjdotr;
            if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) 
              e += -qj*r3inv*pidotr;
            if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0)
              e += qtmp*r3inv*pjdotr;
            e *= fq;
          } else e = (acctyp)0.0;
            e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          } 
        }
        if (vflag>0) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }

    } // for nbor
    store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

__kernel void k_dipole_lj_fast(const __global numtyp4 *restrict x_, 
                               const __global numtyp4 *restrict lj1_in,
                               const __global numtyp4 *restrict lj3_in, 
                               const __global numtyp *restrict sp_lj_in,
                               const __global int *dev_nbor, 
                               const __global int *dev_packed,
                               __global acctyp4 *restrict ans, 
                               __global acctyp *restrict engv, 
                               const int eflag, const int vflag, const int inum, 
                               const int nbor_pitch, 
                               const __global numtyp *restrict q_,
                               const __global numtyp4 *restrict mu_,
                               const __global numtyp *restrict _cutsq, 
                               const numtyp qqrd2e, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    cutsq[tid]=_cutsq[tid];
    if (eflag>0)
      lj3[tid]=lj3_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp4 tor;
  tor.x=(acctyp)0;
  tor.y=(acctyp)0;
  tor.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __syncthreads();
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    numtyp4 mui; fetch4(mui,i,mu_tex); //mu_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      numtyp qj; fetch(qj,j,q_tex);
      numtyp4 muj; fetch4(muj,j,mu_tex); //mu_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, r6inv;
        numtyp rinv, r3inv, r5inv, r7inv;
        numtyp pre1, pre2, pre3, pre4;
        numtyp pdotp, pidotr, pjdotr;
        acctyp4 forcecoul, ticoul;
        acctyp4 force;
        
        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;
        } else force_lj = (numtyp)0.0;

        if (rsq < lj1[mtype].w) {
          rinv = ucl_rsqrt(rsq);

          // charge-charge
          if (qtmp != (numtyp)0.0 && qj != (numtyp)0.0) { 
            r3inv = r2inv*rinv;
            pre1 = qtmp*qj*r3inv;

            forcecoul.x += pre1*delx;
            forcecoul.y += pre1*dely;
            forcecoul.z += pre1*delz;
          }

          // dipole-dipole
          if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0) {
            r3inv = r2inv*rinv; 
            r5inv = r3inv*r2inv;
            r7inv = r5inv*r2inv;
            pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

            pre1 = (numtyp)3.0*r5inv*pdotp - (numtyp)15.0*r7inv*pidotr*pjdotr;
            pre2 = (numtyp)3.0*r5inv*pjdotr;
            pre3 = (numtyp)3.0*r5inv*pidotr;
            pre4 = (numtyp)(-1.0)*r3inv;

            forcecoul.x += pre1*delx + pre2*mui.x + pre3*muj.x;
            forcecoul.y += pre1*dely + pre2*mui.y + pre3*muj.y;
            forcecoul.z += pre1*delz + pre2*mui.z + pre3*muj.z;
	    
            numtyp crossx = pre4 * (mui.y*muj.z - mui.z*muj.y);
            numtyp crossy = pre4 * (mui.z*muj.x - mui.x*muj.z);
            numtyp crossz = pre4 * (mui.x*muj.y - mui.y*muj.x);

            ticoul.x += crossx + pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += crossy + pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += crossz + pre2 * (mui.x*dely - mui.y*delx);
          }

          // dipole-charge
          if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) { 
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pre1 = (numtyp)3.0*qj*r5inv * pidotr;
            pre2 = qj*r3inv;
            
            forcecoul.x += pre2*mui.x - pre1*delx;
            forcecoul.y += pre2*mui.y - pre1*dely;
            forcecoul.z += pre2*mui.z - pre1*delz;
            ticoul.x += pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += pre2 * (mui.x*dely - mui.y*delx);
          }
          
          // charge-dipole
          if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;
            pre1 = (numtyp)3.0*qtmp*r5inv * pjdotr;
            pre2 = qtmp*r3inv;
            
            forcecoul.x += pre1*delx - pre2*muj.x;
            forcecoul.y += pre1*dely - pre2*muj.y;
            forcecoul.z += pre1*delz - pre2*muj.z;
          }
        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
          ticoul.x = ticoul.y = ticoul.z = (acctyp)0;
        }

        numtyp fq = factor_coul*qqrd2e;
        force.x = fq*forcecoul.x + delx*force_lj;
        force.y = fq*forcecoul.y + dely*force_lj;
        force.z = fq*forcecoul.z + delz*force_lj;

        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=fq*ticoul.x;
        tor.y+=fq*ticoul.y;
        tor.z+=fq*ticoul.z;

        if (eflag>0) {
          acctyp e = (acctyp)0;  
          if (rsq < lj1[mtype].w) {
            e = qtmp*qj*rinv;
            if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0)
              e += r3inv*pdotp - (numtyp)3.0*r5inv*pidotr*pjdotr;
            if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) 
              e += -qj*r3inv*pidotr;
            if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0)
              e += qtmp*r3inv*pjdotr;
            e *= fq;
          } else e = (acctyp)0;
          e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          }
        }
        if (vflag>0) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }

    } // for nbor
    store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

