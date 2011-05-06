// **************************************************************************
//                                re_squared.cu
//                             -------------------
//                               W. Michael Brown
//
//  Device code for RE-Squared potential acceleration
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Fri May 06 2011
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifndef RE_SQUARED_CU
#define RE_SQUARED_CU

#ifdef NV_KERNEL
#include "ellipsoid_extra.h"
#endif

#define SBBITS 30
#define NEIGHMASK 0x3FFFFFFF
__inline int sbmask(int j) { return j >> SBBITS & 3; }

__kernel void kernel_ellipsoid(__global numtyp4* x_,__global numtyp4 *q,
                               __global numtyp4* shape, __global numtyp4* well, 
                               __global numtyp *gum, __global numtyp2* sig_eps, 
                               const int ntypes, __global numtyp *lshape, 
                               __global int *dev_nbor, const int stride, 
                               __global acctyp4 *ans, const int astride, 
                               __global acctyp *engv, __global int *err_flag, 
                               const int eflag, const int vflag, const int inum,
                               const int nall, const int t_per_atom) {
/*
  int tid=THREAD_ID_X;
  int ii=mul24((int)BLOCK_ID_X,(int)(BLOCK_SIZE_X)/t_per_atom);
  ii+=tid/t_per_atom;
  int offset=tid%t_per_atom;

  __local numtyp sp_lj[4];
  sp_lj[0]=gum[0];    
  sp_lj[1]=gum[1];    
  sp_lj[2]=gum[2];    
  sp_lj[3]=gum[3];    

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  acctyp4 tor;
  tor.x=(acctyp)0;
  tor.y=(acctyp)0;
  tor.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  if (ii<inum) {
    __global int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=stride;
    int numj=*nbor;
    nbor+=stride;
    __global int *nbor_end=nbor+mul24(stride,numj);
    nbor+=mul24(offset,stride);
    int n_stride=mul24(t_per_atom,stride);
  
    numtyp4 ix=x_[i];
    int itype=ix.w;

    numtyp a1[9];       // Rotation matrix (lab->body)
    numtyp aTe1[9];     // A'*E
    numtyp gamma1[9];   // A'*S^2*A
    numtyp sa1[9];      // S^2*A;
    numtyp lA1_1[9], lA1_2[9], lA1_3[9]; // -A*rotation generator (x,y, or z)
    numtyp lAtwo1_1[9], lAtwo1_2[9], lAtwo1_3[9];  // A'*S^2*lA
    numtyp lAsa1_1[9], lAsa1_2[9], lAsa1_3[9];   // lAtwo+lA'*sa
    numtyp4 ishape=shape[itype]; //MAKE SURE SQUARED
    
    {
      numtyp aTs[9];    // A1'*S1^2
      gpu_quat_to_mat_trans(q,i,a1);
      gpu_transpose_times_diag3(a1,well[itype],aTe1);
      gpu_transpose_times_diag3(a1,ishape,aTs);
      gpu_diag_times3(ishape,a1,sa1);
      gpu_times3(aTs1,a1,gamma1);
      gpu_rotation_generator_x(a1,lA1[0]);
      gpu_rotation_generator_y(a1,lA1[1]);
      gpu_rotation_generator_z(a1,lA1[2]);
      gpu_times3(aTs,lA1_1,lAtwo1_1);
      gpu_transpose_times3(lA1_1,sa1,lAsa1_1);
      gpu_plus3(lAsa1_1,lAtwo1_1,lAsa1_1);
      gpu_times3(aTs,lA1_2,lAtwo1_2);
      gpu_transpose_times3(lA1_2,sa1,lAsa1_2);
      gpu_plus3(lAsa1_2,lAtwo1_2,lAsa1_2);
      gpu_times3(aTs,lA1_3,lAtwo1_3);
      gpu_transpose_times3(lA1_3,sa1,lAsa1_3);
      gpu_plus3(lAsa1_3,lAtwo1_3,lAsa1_3);
    }

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp r12[3];
      r12[0] = jx.x-ix.x;
      r12[1] = jx.y-ix.y;
      r12[2] = jx.z-ix.z;
      numtyp ir = gpu_dot3(r12,r12);
*/
/*
      if (eflag>0)
        energy+=u_r*temp2;
      numtyp temp1 = -eta*u_r*factor_lj;
      if (vflag>0) {
        r12[0]*=-r;
        r12[1]*=-r;
        r12[2]*=-r;
        numtyp ft=temp1*dchi[0]-temp2*dUr[0];
        f.x+=ft;
        virial[0]+=r12[0]*ft;
        ft=temp1*dchi[1]-temp2*dUr[1];
        f.y+=ft;
        virial[1]+=r12[1]*ft;
        virial[3]+=r12[0]*ft;
        ft=temp1*dchi[2]-temp2*dUr[2];
        f.z+=ft;
        virial[2]+=r12[2]*ft;
        virial[4]+=r12[0]*ft;
        virial[5]+=r12[1]*ft;
      } else {
        f.x+=temp1*dchi[0]-temp2*dUr[0];
        f.y+=temp1*dchi[1]-temp2*dUr[1];
        f.z+=temp1*dchi[2]-temp2*dUr[2];
      }

      // Torque on 1
      temp1 = -u_r*eta*factor_lj;
      temp2 = -u_r*chi*factor_lj;
      numtyp temp3 = -chi*eta*factor_lj;
      tor.x+=temp1*tchi[0]+temp2*teta[0]+temp3*tUr[0];
      tor.y+=temp1*tchi[1]+temp2*teta[1]+temp3*tUr[1];
      tor.z+=temp1*tchi[2]+temp2*teta[2]+temp3*tUr[2];
*/ 
/*
    } // for nbor
  } // if ii
  
  // Reduce answers
  if (t_per_atom>1) {
    __local acctyp red_acc[7][BLOCK_PAIR];
    
    red_acc[0][tid]=f.x;
    red_acc[1][tid]=f.y;
    red_acc[2][tid]=f.z;
    red_acc[3][tid]=tor.x;
    red_acc[4][tid]=tor.y;
    red_acc[5][tid]=tor.z;

    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
      if (offset < s) {
        for (int r=0; r<6; r++)
          red_acc[r][tid] += red_acc[r][tid+s];
      }
    }
    
    f.x=red_acc[0][tid];
    f.y=red_acc[1][tid];
    f.z=red_acc[2][tid];
    tor.x=red_acc[3][tid];
    tor.y=red_acc[4][tid];
    tor.z=red_acc[5][tid];

    if (eflag>0 || vflag>0) {
      for (int r=0; r<6; r++)
        red_acc[r][tid]=virial[r];
      red_acc[6][tid]=energy;

      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
        if (offset < s) {
          for (int r=0; r<7; r++)
            red_acc[r][tid] += red_acc[r][tid+s];
        }
      }
    
      for (int r=0; r<6; r++)
        virial[r]=red_acc[r][tid];
      energy=red_acc[6][tid];
    }
  }

  // Store answers
  if (ii<inum && offset==0) {
    __global acctyp *ap1=engv+ii;
    if (eflag>0) {
      *ap1=energy;
      ap1+=astride;
    }
    if (vflag>0) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=astride;
      }
    }
    ans[ii]=f;
    ans[ii+astride]=tor;
  } // if ii
*/
}

#endif

