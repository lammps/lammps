/***************************************************************************
                              gb_gpu_kernel.cu
                             -------------------
                               W. Michael Brown

  Routines that actually perform the force/torque computation

   *** Force Decomposition by Atom Version ***
 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Jun 23 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef GB_GPU_KERNEL
#define GB_GPU_KERNEL

#include "gb_gpu_extra.h"

template <class numtyp>
static __inline__ __device__ void compute_eta_torque(numtyp m[9], 
                                                     numtyp m2[9],
                                                     const int i, 
                                                     numtyp ans[9])
{
  numtyp den = m[3]*m[2]*m[7]-m[0]*m[5]*m[7]-
    m[2]*m[6]*m[4]+m[1]*m[6]*m[5]-
    m[3]*m[1]*m[8]+m[0]*m[4]*m[8];
  den = (numtyp)1.0/den;
  
  numtyp shapex=_shape_<numtyp>(i,0);
  numtyp shapey=_shape_<numtyp>(i,1);
  numtyp shapez=_shape_<numtyp>(i,2);
  
  ans[0] = shapex*(m[5]*m[1]*m2[2]+(numtyp)2.0*m[4]*m[8]*m2[0]-
		    m[4]*m2[2]*m[2]-(numtyp)2.0*m[5]*m2[0]*m[7]+
		    m2[1]*m[2]*m[7]-m2[1]*m[1]*m[8]-
		    m[3]*m[8]*m2[1]+m[6]*m[5]*m2[1]+
		    m[3]*m2[2]*m[7]-m2[2]*m[6]*m[4])*den;
  
  ans[1] = shapex*(m[2]*m2[0]*m[7]-m[8]*m2[0]*m[1]+
		    (numtyp)2.0*m[0]*m[8]*m2[1]-m[0]*m2[2]*m[5]-
		    (numtyp)2.0*m[6]*m[2]*m2[1]+m2[2]*m[3]*m[2]-
		    m[8]*m[3]*m2[0]+m[6]*m2[0]*m[5]+
		    m[6]*m2[2]*m[1]-m2[2]*m[0]*m[7])*den;
  
  ans[2] = shapex*(m[1]*m[5]*m2[0]-m[2]*m2[0]*m[4]-
		    m[0]*m[5]*m2[1]+m[3]*m[2]*m2[1]-
		    m2[1]*m[0]*m[7]-m[6]*m[4]*m2[0]+
		    (numtyp)2.0*m[4]*m[0]*m2[2]-(numtyp)2.0*m[3]*m2[2]*m[1]+
		    m[3]*m[7]*m2[0]+m[6]*m2[1]*m[1])*den;
  
  ans[3] = shapey*(-m[4]*m2[5]*m[2]+(numtyp)2.0*m[4]*m[8]*m2[3]+
		    m[5]*m[1]*m2[5]-(numtyp)2.0*m[5]*m2[3]*m[7]+
		    m2[4]*m[2]*m[7]-m2[4]*m[1]*m[8]-
		    m[3]*m[8]*m2[4]+m[6]*m[5]*m2[4]-
		    m2[5]*m[6]*m[4]+m[3]*m2[5]*m[7])*den;
  
  ans[4] = shapey*(m[2]*m2[3]*m[7]-m[1]*m[8]*m2[3]+
		    (numtyp)2.0*m[8]*m[0]*m2[4]-m2[5]*m[0]*m[5]-
		    (numtyp)2.0*m[6]*m2[4]*m[2]-m[3]*m[8]*m2[3]+
		    m[6]*m[5]*m2[3]+m[3]*m2[5]*m[2]-
		    m[0]*m2[5]*m[7]+m2[5]*m[1]*m[6])*den;
  
  ans[5] = shapey*(m[1]*m[5]*m2[3]-m[2]*m2[3]*m[4]-
		    m[0]*m[5]*m2[4]+m[3]*m[2]*m2[4]+
		    (numtyp)2.0*m[4]*m[0]*m2[5]-m[0]*m2[4]*m[7]+
		    m[1]*m[6]*m2[4]-m2[3]*m[6]*m[4]-
		    (numtyp)2.0*m[3]*m[1]*m2[5]+m[3]*m2[3]*m[7])*den;
  
  ans[6] = shapez*(-m[4]*m[2]*m2[8]+m[1]*m[5]*m2[8]+
		    (numtyp)2.0*m[4]*m2[6]*m[8]-m[1]*m2[7]*m[8]+
		    m[2]*m[7]*m2[7]-(numtyp)2.0*m2[6]*m[7]*m[5]-
		    m[3]*m2[7]*m[8]+m[5]*m[6]*m2[7]-
		    m[4]*m[6]*m2[8]+m[7]*m[3]*m2[8])*den;
  
  ans[7] = shapez*-(m[1]*m[8]*m2[6]-m[2]*m2[6]*m[7]-
		     (numtyp)2.0*m2[7]*m[0]*m[8]+m[5]*m2[8]*m[0]+
		     (numtyp)2.0*m2[7]*m[2]*m[6]+m[3]*m2[6]*m[8]-
		     m[3]*m[2]*m2[8]-m[5]*m[6]*m2[6]+
		     m[0]*m2[8]*m[7]-m2[8]*m[1]*m[6])*den;
  
  ans[8] = shapez*(m[1]*m[5]*m2[6]-m[2]*m2[6]*m[4]-
		    m[0]*m[5]*m2[7]+m[3]*m[2]*m2[7]-
		    m[4]*m[6]*m2[6]-m[7]*m2[7]*m[0]+
		    (numtyp)2.0*m[4]*m2[8]*m[0]+m[7]*m[3]*m2[6]+
		    m[6]*m[1]*m2[7]-(numtyp)2.0*m2[8]*m[3]*m[1])*den;
}

template<class numtyp, class acctyp>
__global__ void kernel_gayberne(const numtyp *gum, const numtyp *special_lj,
                                const int *dev_nbor, const int nbor_pitch, 
                                acctyp *ans, size_t ans_pitch, int *err_flag, 
                                const bool eflag, const bool vflag,
                                const int inum, const int nall) {
                                
  __shared__ numtyp sp_lj[4];

   // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  ii+=INT_MUL(blockIdx.x,blockDim.x);                                  

  if (ii<inum) {
  
  acctyp energy=(numtyp)0;
  acctyp fx=(numtyp)0;
  acctyp fy=(numtyp)0;
  acctyp fz=(numtyp)0;
  acctyp torx=(numtyp)0;
  acctyp tory=(numtyp)0;
  acctyp torz=(numtyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(numtyp)0;
  
  const int *nbor=dev_nbor+ii;
  int i=*nbor;
  nbor+=nbor_pitch;
  nbor+=nbor_pitch;
  nbor+=nbor_pitch;
  int numj=*nbor;
  nbor+=nbor_pitch;
  const int *nbor_end=nbor+INT_MUL(nbor_pitch,numj);
  
  numtyp ix=_x_<numtyp>(i,0);
  numtyp iy=_x_<numtyp>(i,1);
  numtyp iz=_x_<numtyp>(i,2);
  int itype=_x_<numtyp>(i,7);
  numtyp a1[9], b1[9], g1[9];
  {
    numtyp t[9];
    gpu_quat_to_mat_trans(i,a1);
    gpu_shape_times3(itype,a1,t);
    gpu_transpose_times3(a1,t,g1);
    gpu_well_times3(itype,a1,t);
    gpu_transpose_times3(a1,t,b1);
  }

  numtyp factor_lj;
  for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
  
  int j=*nbor;
  if (j < nall) 
    factor_lj = (numtyp)1.0;
  else {
    factor_lj = sp_lj[j/nall];
    j %= nall;
  }
  int jtype=_x_<numtyp>(j,7);

  // Compute r12
  numtyp r12[3];
  r12[0] = _x_<numtyp>(j,0)-ix;
  r12[1] = _x_<numtyp>(j,1)-iy;
  r12[2] = _x_<numtyp>(j,2)-iz;
  numtyp ir = gpu_dot3(r12,r12);

  ir = rsqrt(ir);
  numtyp r = (numtyp)1.0/ir;

  numtyp a2[9];
  gpu_quat_to_mat_trans(j,a2);
  
  numtyp u_r, dUr[3], tUr[3], eta, teta[3];
  { // Compute U_r, dUr, eta, and teta
    // Compute g12
    numtyp g12[9];
    {
      numtyp g2[9];
      {
          gpu_shape_times3(jtype,a2,g12);
          gpu_transpose_times3(a2,g12,g2);
          gpu_plus3(g1,g2,g12);
      }
  
      { // Compute U_r and dUr
    
        // Compute kappa
        numtyp kappa[3];
        gpu_mldivide3(g12,r12,kappa,err_flag);

        // -- replace r12 with r12 hat
        r12[0]*=ir;
        r12[1]*=ir;
        r12[2]*=ir;

        // -- kappa is now / r
        kappa[0]*=ir;
        kappa[1]*=ir;
        kappa[2]*=ir;
  
        // energy
  
        // compute u_r and dUr
        numtyp uslj_rsq;
        {
          // Compute distance of closest approach
          numtyp h12, sigma12;
          sigma12 = gpu_dot3(r12,kappa);
          sigma12 = rsqrt((numtyp)0.5*sigma12);
          h12 = r-sigma12;

          // -- kappa is now ok
          kappa[0]*=r;
          kappa[1]*=r;
          kappa[2]*=r;
          
          numtyp sigma = _sigma_<numtyp>(itype,jtype);
          numtyp epsilon = _epsilon_<numtyp>(itype,jtype);
          numtyp varrho = sigma/(h12+gum[0]*sigma);
          numtyp varrho6 = varrho*varrho*varrho;
          varrho6*=varrho6;
          numtyp varrho12 = varrho6*varrho6;
          u_r = (numtyp)4.0*epsilon*(varrho12-varrho6);

          numtyp temp1 = ((numtyp)2.0*varrho12*varrho-varrho6*varrho)/sigma;
          temp1 = temp1*(numtyp)24.0*epsilon;
          uslj_rsq = temp1*sigma12*sigma12*sigma12*(numtyp)0.5;
          numtyp temp2 = gpu_dot3(kappa,r12);
          uslj_rsq = uslj_rsq*ir*ir;

          dUr[0] = temp1*r12[0]+uslj_rsq*(kappa[0]-temp2*r12[0]);
          dUr[1] = temp1*r12[1]+uslj_rsq*(kappa[1]-temp2*r12[1]);
          dUr[2] = temp1*r12[2]+uslj_rsq*(kappa[2]-temp2*r12[2]);
        }

        // torque for particle 1
        {
          numtyp tempv[3], tempv2[3];
          tempv[0] = -uslj_rsq*kappa[0];
          tempv[1] = -uslj_rsq*kappa[1];
          tempv[2] = -uslj_rsq*kappa[2];
          gpu_row_times3(kappa,g1,tempv2);
          gpu_cross3(tempv,tempv2,tUr);
        }
      }
    }
     
    // Compute eta
    {
      eta = (numtyp)2.0*_lshape_<numtyp>(itype)*_lshape_<numtyp>(jtype);
      numtyp det_g12 = gpu_det3(g12);
      eta = pow(eta/det_g12,gum[1]);
    }
    
    // Compute teta
    numtyp temp[9], tempv[3], tempv2[3];
    compute_eta_torque(g12,a1,itype,temp);
    numtyp temp1 = -eta*gum[1];

    tempv[0] = temp1*temp[0];
    tempv[1] = temp1*temp[1];
    tempv[2] = temp1*temp[2];
    gpu_cross3(a1,tempv,tempv2);
    teta[0] = tempv2[0];
    teta[1] = tempv2[1];
    teta[2] = tempv2[2];
  
    tempv[0] = temp1*temp[3];
    tempv[1] = temp1*temp[4];
    tempv[2] = temp1*temp[5];
    gpu_cross3(a1+3,tempv,tempv2);
    teta[0] += tempv2[0];
    teta[1] += tempv2[1];
    teta[2] += tempv2[2];

    tempv[0] = temp1*temp[6];
    tempv[1] = temp1*temp[7];
    tempv[2] = temp1*temp[8];
    gpu_cross3(a1+6,tempv,tempv2);
    teta[0] += tempv2[0];
    teta[1] += tempv2[1];
    teta[2] += tempv2[2];
  }
  
  numtyp chi, dchi[3], tchi[3];
  { // Compute chi and dchi

    // Compute b12
    numtyp b2[9], b12[9];
    {
      gpu_well_times3(jtype,a2,b12);
      gpu_transpose_times3(a2,b12,b2);
      gpu_plus3(b1,b2,b12);
    }

    // compute chi_12
    r12[0]*=r;
    r12[1]*=r;
    r12[2]*=r;
    numtyp iota[3];
    gpu_mldivide3(b12,r12,iota,err_flag);
    // -- iota is now iota/r
    iota[0]*=ir;
    iota[1]*=ir;
    iota[2]*=ir;
    r12[0]*=ir;
    r12[1]*=ir;
    r12[2]*=ir;
    chi = gpu_dot3(r12,iota);
    chi = pow(chi*(numtyp)2.0,gum[2]);

    // -- iota is now ok
    iota[0]*=r;
    iota[1]*=r;
    iota[2]*=r;

    numtyp temp1 = gpu_dot3(iota,r12);
    numtyp temp2 = (numtyp)-4.0*ir*ir*gum[2]*pow(chi,(gum[2]-(numtyp)1.0)/gum[2]);
    dchi[0] = temp2*(iota[0]-temp1*r12[0]);
    dchi[1] = temp2*(iota[1]-temp1*r12[1]);
    dchi[2] = temp2*(iota[2]-temp1*r12[2]);

    // compute t_chi
    numtyp tempv[3];
    gpu_row_times3(iota,b1,tempv);
    gpu_cross3(tempv,iota,tchi);
    temp1 = (numtyp)-4.0*ir*ir;
    tchi[0] *= temp1;
    tchi[1] *= temp1;
    tchi[2] *= temp1;
  }

  numtyp temp2 = factor_lj*eta*chi;
  if (eflag)
    energy+=u_r*temp2;
  numtyp temp1 = -eta*u_r*factor_lj;
  if (vflag) {
    r12[0]*=-r;
    r12[1]*=-r;
    r12[2]*=-r;
    numtyp ft=temp1*dchi[0]-temp2*dUr[0];
    fx+=ft;
    virial[0]+=r12[0]*ft;
    ft=temp1*dchi[1]-temp2*dUr[1];
    fy+=ft;
    virial[1]+=r12[1]*ft;
    virial[3]+=r12[0]*ft;
    ft=temp1*dchi[2]-temp2*dUr[2];
    fz+=ft;
    virial[2]+=r12[2]*ft;
    virial[4]+=r12[0]*ft;
    virial[5]+=r12[1]*ft;
  } else {
    fx+=temp1*dchi[0]-temp2*dUr[0];
    fy+=temp1*dchi[1]-temp2*dUr[1];
    fz+=temp1*dchi[2]-temp2*dUr[2];
  }

  // Torque on 1
  temp1 = -u_r*eta*factor_lj;
  temp2 = -u_r*chi*factor_lj;
  numtyp temp3 = -chi*eta*factor_lj;
  torx+=temp1*tchi[0]+temp2*teta[0]+temp3*tUr[0];
  tory+=temp1*tchi[1]+temp2*teta[1]+temp3*tUr[1];
  torz+=temp1*tchi[2]+temp2*teta[2]+temp3*tUr[2];

  } // for nbor

  // Store answers
  acctyp *ap1=ans+ii;
  if (eflag) {
    *ap1=energy;
    ap1+=ans_pitch;
  }
  if (vflag) {
    for (int i=0; i<6; i++) {
      *ap1=virial[i];
      ap1+=ans_pitch;
    }
  }
  *ap1=fx;
  ap1+=ans_pitch;
  *ap1=fy;
  ap1+=ans_pitch;
  *ap1=fz;
  ap1+=ans_pitch;
  *ap1=torx;
  ap1+=ans_pitch;
  *ap1=tory;
  ap1+=ans_pitch;
  *ap1=torz;

  } // if ii
}

template<class numtyp, class acctyp>
__global__ void kernel_sphere_gb(const numtyp *gum, const numtyp *special_lj,
                                 const int *dev_nbor, const int nbor_pitch, 
                                 acctyp *ans, size_t ans_pitch, int *err_flag, 
                                 const bool eflag, const bool vflag,
                                 const int start, const int inum, 
                                 const int nall) {
  __shared__ numtyp sp_lj[4];

  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  ii+=INT_MUL(blockIdx.x,blockDim.x)+start;

  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    nbor+=nbor_pitch;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *nbor_end=nbor+INT_MUL(nbor_pitch,numj);
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=_x_<numtyp>(i,7);
      
    numtyp oner=_shape_<numtyp>(itype,0);
    numtyp one_well=_well_<numtyp>(itype,0);
  
    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
  
      int j=*nbor;
      if (j < nall) 
        factor_lj = (numtyp)1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int jtype=_x_<numtyp>(j,7);

      // Compute r12
      numtyp r12[3];
      r12[0] = _x_<numtyp>(j,0)-ix;
      r12[1] = _x_<numtyp>(j,1)-iy;
      r12[2] = _x_<numtyp>(j,2)-iz;
      numtyp ir = gpu_dot3(r12,r12);

      ir = rsqrt(ir);
      numtyp r = (numtyp)1.0/ir;
      
      numtyp r12hat[3];
      r12hat[0]=r12[0]*ir;
      r12hat[1]=r12[1]*ir;
      r12hat[2]=r12[2]*ir;

      numtyp a2[9];
      gpu_quat_to_mat_trans(j,a2);
  
      numtyp u_r, dUr[3], eta;
      { // Compute U_r, dUr, eta, and teta
        // Compute g12
        numtyp g12[9];
        {
          {
            numtyp g2[9];
            gpu_shape_times3(jtype,a2,g12);
            gpu_transpose_times3(a2,g12,g2);
            g12[0]=g2[0]+oner;
            g12[4]=g2[4]+oner;
            g12[8]=g2[8]+oner;
            g12[1]=g2[1];
            g12[2]=g2[2];
            g12[3]=g2[3];
            g12[5]=g2[5];
            g12[6]=g2[6];
            g12[7]=g2[7];    
          }
  
          { // Compute U_r and dUr
    
            // Compute kappa
            numtyp kappa[3];
            gpu_mldivide3(g12,r12,kappa,err_flag);

            // -- kappa is now / r
            kappa[0]*=ir;
            kappa[1]*=ir;
            kappa[2]*=ir;
  
            // energy
  
            // compute u_r and dUr
            numtyp uslj_rsq;
            {
              // Compute distance of closest approach
              numtyp h12, sigma12;
              sigma12 = gpu_dot3(r12hat,kappa);
              sigma12 = rsqrt((numtyp)0.5*sigma12);
              h12 = r-sigma12;

              // -- kappa is now ok
              kappa[0]*=r;
              kappa[1]*=r;
              kappa[2]*=r;
          
              numtyp sigma = _sigma_<numtyp>(itype,jtype);
              numtyp epsilon = _epsilon_<numtyp>(itype,jtype);
              numtyp varrho = sigma/(h12+gum[0]*sigma);
              numtyp varrho6 = varrho*varrho*varrho;
              varrho6*=varrho6;
              numtyp varrho12 = varrho6*varrho6;
              u_r = (numtyp)4.0*epsilon*(varrho12-varrho6);

              numtyp temp1 = ((numtyp)2.0*varrho12*varrho-varrho6*varrho)/sigma;
              temp1 = temp1*(numtyp)24.0*epsilon;
              uslj_rsq = temp1*sigma12*sigma12*sigma12*(numtyp)0.5;
              numtyp temp2 = gpu_dot3(kappa,r12hat);
              uslj_rsq = uslj_rsq*ir*ir;

              dUr[0] = temp1*r12hat[0]+uslj_rsq*(kappa[0]-temp2*r12hat[0]);
              dUr[1] = temp1*r12hat[1]+uslj_rsq*(kappa[1]-temp2*r12hat[1]);
              dUr[2] = temp1*r12hat[2]+uslj_rsq*(kappa[2]-temp2*r12hat[2]);
            }
          }
        }
     
        // Compute eta
        {
          eta = (numtyp)2.0*_lshape_<numtyp>(itype)*_lshape_<numtyp>(jtype);
          numtyp det_g12 = gpu_det3(g12);
          eta = pow(eta/det_g12,gum[1]);
        }
      }
  
      numtyp chi, dchi[3];
      { // Compute chi and dchi

        // Compute b12
        numtyp b12[9];
        {
          numtyp b2[9];
          gpu_well_times3(jtype,a2,b12);
          gpu_transpose_times3(a2,b12,b2);
          b12[0]=b2[0]+one_well;
          b12[4]=b2[4]+one_well;
          b12[8]=b2[8]+one_well;
          b12[1]=b2[1];
          b12[2]=b2[2];
          b12[3]=b2[3];
          b12[5]=b2[5];
          b12[6]=b2[6];
          b12[7]=b2[7];    
        }

        // compute chi_12
        numtyp iota[3];
        gpu_mldivide3(b12,r12,iota,err_flag);
        // -- iota is now iota/r
        iota[0]*=ir;
        iota[1]*=ir;
        iota[2]*=ir;
        chi = gpu_dot3(r12hat,iota);
        chi = pow(chi*(numtyp)2.0,gum[2]);

        // -- iota is now ok
        iota[0]*=r;
        iota[1]*=r;
        iota[2]*=r;

        numtyp temp1 = gpu_dot3(iota,r12hat);
        numtyp temp2 = (numtyp)-4.0*ir*ir*gum[2]*pow(chi,(gum[2]-(numtyp)1.0)/gum[2]);
        dchi[0] = temp2*(iota[0]-temp1*r12hat[0]);
        dchi[1] = temp2*(iota[1]-temp1*r12hat[1]);
        dchi[2] = temp2*(iota[2]-temp1*r12hat[2]);
      }

      numtyp temp2 = factor_lj*eta*chi;
      if (eflag)
        energy+=u_r*temp2;
      numtyp temp1 = -eta*u_r*factor_lj;
      if (vflag) {
        r12[0]*=-1;
        r12[1]*=-1;
        r12[2]*=-1;
        numtyp ft=temp1*dchi[0]-temp2*dUr[0];
        fx+=ft;
        virial[0]+=r12[0]*ft;
        ft=temp1*dchi[1]-temp2*dUr[1];
        fy+=ft;
        virial[1]+=r12[1]*ft;
        virial[3]+=r12[0]*ft;
        ft=temp1*dchi[2]-temp2*dUr[2];
        fz+=ft;
        virial[2]+=r12[2]*ft;
        virial[4]+=r12[0]*ft;
        virial[5]+=r12[1]*ft;
      } else {
        fx+=temp1*dchi[0]-temp2*dUr[0];
        fy+=temp1*dchi[1]-temp2*dUr[1];
        fz+=temp1*dchi[2]-temp2*dUr[2];
      }
    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1=fx;
    ap1+=ans_pitch;
    *ap1=fy;
    ap1+=ans_pitch;
    *ap1=fz;
  } // if ii
}

template<class numtyp, class acctyp>
__global__ void kernel_lj(const numtyp *special_lj, const int *dev_nbor, 
                          const int *dev_ij, const int nbor_pitch, acctyp *ans, 
                          size_t ans_pitch, int *err_flag, const bool eflag, 
                          const bool vflag, const int start, const int inum,
                          const int nall) {
  __shared__ numtyp sp_lj[4];                              
  
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  ii+=INT_MUL(blockIdx.x,blockDim.x)+start;

  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=_x_<numtyp>(i,7);

    numtyp factor_lj;
    for ( ; list<list_end; list++) {
  
      int j=*list;
      if (j < nall) 
        factor_lj = (numtyp)1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int jtype=_x_<numtyp>(j,7);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      if (r2inv<_cutsq_<numtyp>(itype,jtype) &&
          _form_(itype,jtype)==SPHERE_SPHERE) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = r2inv*r6inv*(_lj1_<numtyp>(itype,jtype).x*r6inv-
                                    _lj1_<numtyp>(itype,jtype).y);
        force*=factor_lj;
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(_lj3_<numtyp>(itype,jtype).x*r6inv-
                          _lj3_<numtyp>(itype,jtype).y);
          energy+=factor_lj*(e-_offset_<numtyp>(1,1)); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1+=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1+=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1+=fx;
    ap1+=ans_pitch;
    *ap1+=fy;
    ap1+=ans_pitch;
    *ap1+=fz;

  } // if ii
}

template<class numtyp, class acctyp>
__global__ void kernel_lj_fast(const numtyp *special_lj, const int *dev_nbor, 
                               const int *dev_ij, const int nbor_pitch, 
                               acctyp *ans, size_t ans_pitch,int *err_flag,
                               const bool eflag, const bool vflag, 
                               const int start, const int inum, const int nall){
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  __shared__ numtyp sp_lj[4];                              
  __shared__ int form[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj4[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp offset[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    int itype=ii/MAX_SHARED_TYPES;
    int jtype=ii%MAX_SHARED_TYPES;
    cutsq[ii]=_cutsq_<numtyp>(itype,jtype);
    form[ii]=_form_(itype,jtype);
    lj1[ii]=_lj1_<numtyp>(itype,jtype).x;
    lj2[ii]=_lj1_<numtyp>(itype,jtype).y;
    if (eflag) {
      lj3[ii]=_lj3_<numtyp>(itype,jtype).x;
      lj4[ii]=_lj3_<numtyp>(itype,jtype).y;
      offset[ii]=_offset_<numtyp>(itype,jtype);
    }
  }
  ii+=INT_MUL(blockIdx.x,blockDim.x)+start;
  
  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=INT_MUL(MAX_SHARED_TYPES,_x_<numtyp>(i,7));

    numtyp factor_lj;
    for ( ; list<list_end; list++) {
  
      int j=*list;
      if (j < nall) 
        factor_lj = (numtyp)1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int mtype=itype+_x_<numtyp>(j,7);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      if (r2inv<cutsq[mtype] && form[mtype]==SPHERE_SPHERE) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = factor_lj*r2inv*r6inv*(lj1[mtype]*r6inv-lj2[mtype]);
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(lj3[mtype]*r6inv-lj4[mtype]);
          energy+=factor_lj*(e-offset[mtype]); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1+=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1+=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1+=fx;
    ap1+=ans_pitch;
    *ap1+=fy;
    ap1+=ans_pitch;
    *ap1+=fz;

  } // if ii
}

#endif
