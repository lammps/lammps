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
      numtyp r[3], rhat[3];
      r[0] = jx.x-ix.x;
      r[1] = jx.y-ix.y;
      r[2] = jx.z-ix.z;
      numtyp rnorm = gpu_dot3(r12,r12);
      rnorm = sqrt(rnorm);
      rhat[0] = r[0]/rnorm;
      rhat[1] = r[1]/rnorm;
      rhat[2] = r[2]/rnorm;

      numtyp temp[9], s[3], z1[3], z2[3], v1[3], v2[3];
      numtyp sigma12, sigma1, sigma2;
      gpu_plus3(gamma1,gamma2,temp);
      gpu_mldivide3(temp,rhat,s,error);
      sigma12 = 1.0/sqrt(0.5*gpu_dot3(s,rhat));
      gpu_times_column3(A1,rhat,z1);
      gpu_times_column3(A2,rhat,z2);
      v1[0] = z1[0]/shape2[type[i]].x;
      v1[1] = z1[1]/shape2[type[i]].y;
      v1[2] = z1[2]/shape2[type[i]].z;
      v2[0] = z2[0]/shape2[type[j]].x;
      v2[1] = z2[1]/shape2[type[j]].y;
      v2[2] = z2[2]/shape2[type[j]].z;
      sigma1 = 1.0/sqrt(gpu_dot3(z1,v1));
      sigma2 = 1.0/sqrt(gpu_dot3(z2,v2));
  H12[0][0] = wi.gamma[0][0]/sigma1+wj.gamma[0][0]/sigma2;
  H12[0][1] = wi.gamma[0][1]/sigma1+wj.gamma[0][1]/sigma2;
  H12[0][2] = wi.gamma[0][2]/sigma1+wj.gamma[0][2]/sigma2;
  H12[1][0] = wi.gamma[1][0]/sigma1+wj.gamma[1][0]/sigma2;
  H12[1][1] = wi.gamma[1][1]/sigma1+wj.gamma[1][1]/sigma2;
  H12[1][2] = wi.gamma[1][2]/sigma1+wj.gamma[1][2]/sigma2;
  H12[2][0] = wi.gamma[2][0]/sigma1+wj.gamma[2][0]/sigma2;
  H12[2][1] = wi.gamma[2][1]/sigma1+wj.gamma[2][1]/sigma2;
  H12[2][2] = wi.gamma[2][2]/sigma1+wj.gamma[2][2]/sigma2;
  dH=MathExtra::det3(H12);
  sigma1p2 = sigma1*sigma1;
  sigma2p2 = sigma2*sigma2;
  lambda = lshape[type[i]]/sigma1p2 + lshape[type[j]]/sigma2p2;
  nu = sqrt(dH/(sigma1+sigma2));
  MathExtra::times3(wi.aTe,wi.A,temp);
  double temp2[3][3];
  MathExtra::times3(wj.aTe,wj.A,temp2);
  MathExtra::plus3(temp,temp2,temp);
  MathExtra::mldivide3(temp,rhat,w,error);
  h12 = rnorm-sigma12;
  eta = lambda/nu;
  chi = 2.0*MathExtra::dot3(rhat,w);
  sprod = lshape[type[i]] * lshape[type[j]];
  sigh = sigma[type[i]][type[j]]/h12;
  tprod = eta*chi*sigh;

  double stemp = h12/2.0;
  Ua = (shape1[type[i]][0]+stemp)*(shape1[type[i]][1]+stemp)*
       (shape1[type[i]][2]+stemp)*(shape1[type[j]][0]+stemp)*
       (shape1[type[j]][1]+stemp)*(shape1[type[j]][2]+stemp);
  Ua = (1.0+3.0*tprod)*sprod/Ua;
  Ua = epsilon[type[i]][type[j]]*Ua/-36.0;

  stemp = h12/cr60;
  Ur = (shape1[type[i]][0]+stemp)*(shape1[type[i]][1]+stemp)*
       (shape1[type[i]][2]+stemp)*(shape1[type[j]][0]+stemp)*
       (shape1[type[j]][1]+stemp)*(shape1[type[j]][2]+stemp);
  Ur = (1.0+b_alpha*tprod)*sprod/Ur;
  Ur = epsilon[type[i]][type[j]]*Ur*pow(sigh,6.0)/2025.0;

  // force

  sec = sigma[type[i]][type[j]]*eta*chi;
  sigma12p3 = pow(sigma12,3.0);
  sigma1p3 = sigma1p2*sigma1;
  sigma2p3 = sigma2p2*sigma2;
  vsigma1[0] = -sigma1p3*v1[0];
  vsigma1[1] = -sigma1p3*v1[1];
  vsigma1[2] = -sigma1p3*v1[2];
  vsigma2[0] = -sigma2p3*v2[0];
  vsigma2[1] = -sigma2p3*v2[1];
  vsigma2[2] = -sigma2p3*v2[2];
  gsigma1[0][0] = -wi.gamma[0][0]/sigma1p2;
  gsigma1[0][1] = -wi.gamma[0][1]/sigma1p2;
  gsigma1[0][2] = -wi.gamma[0][2]/sigma1p2;
  gsigma1[1][0] = -wi.gamma[1][0]/sigma1p2;
  gsigma1[1][1] = -wi.gamma[1][1]/sigma1p2;
  gsigma1[1][2] = -wi.gamma[1][2]/sigma1p2;
  gsigma1[2][0] = -wi.gamma[2][0]/sigma1p2;
  gsigma1[2][1] = -wi.gamma[2][1]/sigma1p2;
  gsigma1[2][2] = -wi.gamma[2][2]/sigma1p2;
  gsigma2[0][0] = -wj.gamma[0][0]/sigma2p2;
  gsigma2[0][1] = -wj.gamma[0][1]/sigma2p2;
  gsigma2[0][2] = -wj.gamma[0][2]/sigma2p2;
  gsigma2[1][0] = -wj.gamma[1][0]/sigma2p2;
  gsigma2[1][1] = -wj.gamma[1][1]/sigma2p2;
  gsigma2[1][2] = -wj.gamma[1][2]/sigma2p2;
  gsigma2[2][0] = -wj.gamma[2][0]/sigma2p2;
  gsigma2[2][1] = -wj.gamma[2][1]/sigma2p2;
  gsigma2[2][2] = -wj.gamma[2][2]/sigma2p2;
  tsig1sig2 = eta/(2.0*(sigma1+sigma2));
  tdH = eta/(2.0*dH);
  teta1 = 2.0*eta/lambda;
  teta2 = teta1*lshape[type[j]]/sigma2p3;
  teta1 = teta1*lshape[type[i]]/sigma1p3;
  fourw[0] = 4.0*w[0];
  fourw[1] = 4.0*w[1];
  fourw[2] = 4.0*w[2];
  spr[0] = 0.5*sigma12p3*s[0];
  spr[1] = 0.5*sigma12p3*s[1];
  spr[2] = 0.5*sigma12p3*s[2];

  stemp = 1.0/(shape1[type[i]][0]*2.0+h12)+
          1.0/(shape1[type[i]][1]*2.0+h12)+
          1.0/(shape1[type[i]][2]*2.0+h12)+
          1.0/(shape1[type[j]][0]*2.0+h12)+
          1.0/(shape1[type[j]][1]*2.0+h12)+
          1.0/(shape1[type[j]][2]*2.0+h12);
  hsec = h12+3.0*sec;
  dspu = 1.0/h12-1.0/hsec+stemp;
  pbsu = 3.0*sigma[type[i]][type[j]]/hsec;
  
  stemp = 1.0/(shape1[type[i]][0]*cr60+h12)+
          1.0/(shape1[type[i]][1]*cr60+h12)+
          1.0/(shape1[type[i]][2]*cr60+h12)+
          1.0/(shape1[type[j]][0]*cr60+h12)+
          1.0/(shape1[type[j]][1]*cr60+h12)+
          1.0/(shape1[type[j]][2]*cr60+h12);
  hsec = h12+b_alpha*sec;
  dspr = 7.0/h12-1.0/hsec+stemp;
  pbsr = b_alpha*sigma[type[i]][type[j]]/hsec;
  
  for (int i=0; i<3; i++) {
    u[0] = -rhat[i]*rhat[0];
    u[1] = -rhat[i]*rhat[1];
    u[2] = -rhat[i]*rhat[2];
    u[i] += 1.0;
    u[0] /= rnorm;
    u[1] /= rnorm;
    u[2] /= rnorm;
    MathExtra::times_column3(wi.A,u,u1);
    MathExtra::times_column3(wj.A,u,u2);
    dsigma1=MathExtra::dot3(u1,vsigma1);
    dsigma2=MathExtra::dot3(u2,vsigma2);
    dH12[0][0] = dsigma1*gsigma1[0][0]+dsigma2*gsigma2[0][0];
    dH12[0][1] = dsigma1*gsigma1[0][1]+dsigma2*gsigma2[0][1];
    dH12[0][2] = dsigma1*gsigma1[0][2]+dsigma2*gsigma2[0][2];
    dH12[1][0] = dsigma1*gsigma1[1][0]+dsigma2*gsigma2[1][0];
    dH12[1][1] = dsigma1*gsigma1[1][1]+dsigma2*gsigma2[1][1];
    dH12[1][2] = dsigma1*gsigma1[1][2]+dsigma2*gsigma2[1][2];
    dH12[2][0] = dsigma1*gsigma1[2][0]+dsigma2*gsigma2[2][0];
    dH12[2][1] = dsigma1*gsigma1[2][1]+dsigma2*gsigma2[2][1];
    dH12[2][2] = dsigma1*gsigma1[2][2]+dsigma2*gsigma2[2][2];
    ddH = det_prime(H12,dH12);
    deta = (dsigma1+dsigma2)*tsig1sig2;
    deta -= ddH*tdH;
    deta -= dsigma1*teta1+dsigma2*teta2;
    dchi = MathExtra::dot3(u,fourw);
    dh12 = rhat[i]+MathExtra::dot3(u,spr);
    dUa = pbsu*(eta*dchi+deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi+deta*chi)-dh12*dspr;
    fforce[i]=dUr*Ur+dUa*Ua;
  }
    
  // torque on i

  MathExtra::row_times3(fourw,wi.aTe,fwae);

  for (int i=0; i<3; i++) {
    MathExtra::times_column3(wi.lA[i],rhat,p);
    dsigma1 = MathExtra::dot3(p,vsigma1);
    dH12[0][0] = wi.lAsa[i][0][0]/sigma1+dsigma1*gsigma1[0][0];
    dH12[0][1] = wi.lAsa[i][0][1]/sigma1+dsigma1*gsigma1[0][1];
    dH12[0][2] = wi.lAsa[i][0][2]/sigma1+dsigma1*gsigma1[0][2];
    dH12[1][0] = wi.lAsa[i][1][0]/sigma1+dsigma1*gsigma1[1][0];
    dH12[1][1] = wi.lAsa[i][1][1]/sigma1+dsigma1*gsigma1[1][1];
    dH12[1][2] = wi.lAsa[i][1][2]/sigma1+dsigma1*gsigma1[1][2];
    dH12[2][0] = wi.lAsa[i][2][0]/sigma1+dsigma1*gsigma1[2][0];
    dH12[2][1] = wi.lAsa[i][2][1]/sigma1+dsigma1*gsigma1[2][1];
    dH12[2][2] = wi.lAsa[i][2][2]/sigma1+dsigma1*gsigma1[2][2];
    ddH = det_prime(H12,dH12);
    deta = tsig1sig2*dsigma1-tdH*ddH;
    deta -= teta1*dsigma1;
    double tempv[3];
    MathExtra::times_column3(wi.lA[i],w,tempv);
    dchi = -MathExtra::dot3(fwae,tempv);
    MathExtra::times_column3(wi.lAtwo[i],spr,tempv);
    dh12 = -MathExtra::dot3(s,tempv);

    dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
    ttor[i] = -(dUa*Ua+dUr*Ur);
  }
  
  // torque on j

  if (!(force->newton_pair || j < atom->nlocal))
    return Ua+Ur;

  MathExtra::row_times3(fourw,wj.aTe,fwae);

  for (int i=0; i<3; i++) {
    MathExtra::times_column3(wj.lA[i],rhat,p);
    dsigma2 = MathExtra::dot3(p,vsigma2);
    dH12[0][0] = wj.lAsa[i][0][0]/sigma2+dsigma2*gsigma2[0][0];
    dH12[0][1] = wj.lAsa[i][0][1]/sigma2+dsigma2*gsigma2[0][1];
    dH12[0][2] = wj.lAsa[i][0][2]/sigma2+dsigma2*gsigma2[0][2];
    dH12[1][0] = wj.lAsa[i][1][0]/sigma2+dsigma2*gsigma2[1][0];
    dH12[1][1] = wj.lAsa[i][1][1]/sigma2+dsigma2*gsigma2[1][1];
    dH12[1][2] = wj.lAsa[i][1][2]/sigma2+dsigma2*gsigma2[1][2];
    dH12[2][0] = wj.lAsa[i][2][0]/sigma2+dsigma2*gsigma2[2][0];
    dH12[2][1] = wj.lAsa[i][2][1]/sigma2+dsigma2*gsigma2[2][1];
    dH12[2][2] = wj.lAsa[i][2][2]/sigma2+dsigma2*gsigma2[2][2];
    ddH = det_prime(H12,dH12);
    deta = tsig1sig2*dsigma2-tdH*ddH;
    deta -= teta2*dsigma2;
    double tempv[3];
    MathExtra::times_column3(wj.lA[i],w,tempv);
    dchi = -MathExtra::dot3(fwae,tempv);
    MathExtra::times_column3(wj.lAtwo[i],spr,tempv);
    dh12 = -MathExtra::dot3(s,tempv);

    dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
    rtor[i] = -(dUa*Ua+dUr*Ur);
  }

  return Ua+Ur;
}

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

}

#endif

