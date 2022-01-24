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
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_ellipsoid_extra.h"
#endif

ucl_inline numtyp det_prime(const numtyp m[9], const numtyp m2[9])
{
  numtyp ans;
  ans = m2[0]*m[4]*m[8] - m2[0]*m[5]*m[7] -
        m[3]*m2[1]*m[8] + m[3]*m2[2]*m[7] +
        m[6]*m2[1]*m[5] - m[6]*m2[2]*m[4] +
        m[0]*m2[4]*m[8] - m[0]*m2[5]*m[7] -
        m2[3]*m[1]*m[8] + m2[3]*m[2]*m[7] +
        m[6]*m[1]*m2[5] - m[6]*m[2]*m2[4] +
        m[0]*m[4]*m2[8] - m[0]*m[5]*m2[7] -
        m[3]*m[1]*m2[8] + m[3]*m[2]*m2[7] +
        m2[6]*m[1]*m[5] - m2[6]*m[2]*m[4];
  return ans;
}

__kernel void k_resquared(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict q,
                          const __global numtyp4 *restrict shape,
                          const __global numtyp4 *restrict well,
                          const __global numtyp *restrict splj,
                          const __global numtyp2 *restrict sig_eps,
                          const int ntypes,
                          const __global int *dev_nbor,
                          const int stride,
                          __global acctyp4 *restrict ans,
                          const int astride,
                          __global acctyp *restrict engv,
                          __global int *restrict err_flag,
                          const int eflag, const int vflag, const int inum,
                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[4];
  int n_stride;
  local_allocate_store_ellipse();

  sp_lj[0]=splj[0];
  sp_lj[1]=splj[1];
  sp_lj[2]=splj[2];
  sp_lj[3]=splj[3];

  const numtyp b_alpha=(numtyp)45.0/(numtyp)56.0;
  const numtyp cr60=ucl_cbrt((numtyp)60.0);

  acctyp4 f, tor;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  tor.x=(acctyp)0; tor.y=(acctyp)0; tor.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info_p(dev_nbor,stride,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex);
    int itype=ix.w;

    numtyp a1[9];       // Rotation matrix (lab->body)
    numtyp aTe1[9];     // A'*E
    numtyp gamma1[9];   // A'*S^2*A
    numtyp sa1[9];      // S^2*A;
    numtyp lA1_0[9], lA1_1[9], lA1_2[9]; // -A*rotation generator (x,y, or z)
    numtyp lAtwo1_0[9], lAtwo1_1[9], lAtwo1_2[9];  // A'*S^2*lA
    numtyp lAsa1_0[9], lAsa1_1[9], lAsa1_2[9];   // lAtwo+lA'*sa
    numtyp4 ishape;

    ishape=shape[itype];
    numtyp4 ishape2;
    ishape2.x=ishape.x*ishape.x;
    ishape2.y=ishape.y*ishape.y;
    ishape2.z=ishape.z*ishape.z;
    numtyp ilshape = ishape.x*ishape.y*ishape.z;

    {
      numtyp aTs[9];    // A1'*S1^2
      gpu_quat_to_mat_trans(q,i,a1);
      gpu_transpose_times_diag3(a1,well[itype],aTe1);
      gpu_transpose_times_diag3(a1,ishape2,aTs);
      gpu_diag_times3(ishape2,a1,sa1);
      gpu_times3(aTs,a1,gamma1);
      gpu_rotation_generator_x(a1,lA1_0);
      gpu_rotation_generator_y(a1,lA1_1);
      gpu_rotation_generator_z(a1,lA1_2);
      gpu_times3(aTs,lA1_0,lAtwo1_0);
      gpu_transpose_times3(lA1_0,sa1,lAsa1_0);
      gpu_plus3(lAsa1_0,lAtwo1_0,lAsa1_0);
      gpu_times3(aTs,lA1_1,lAtwo1_1);
      gpu_transpose_times3(lA1_1,sa1,lAsa1_1);
      gpu_plus3(lAsa1_1,lAtwo1_1,lAsa1_1);
      gpu_times3(aTs,lA1_2,lAtwo1_2);
      gpu_transpose_times3(lA1_2,sa1,lAsa1_2);
      gpu_plus3(lAsa1_2,lAtwo1_2,lAsa1_2);
    }
    ishape2.x=ucl_recip(ishape2.x);
    ishape2.y=ucl_recip(ishape2.y);
    ishape2.z=ucl_recip(ishape2.z);

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_nbor[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex);
      int jtype=jx.w;

      // Compute r12
      numtyp r[3], rhat[3];
      numtyp rnorm;
      r[0] = jx.x-ix.x;
      r[1] = jx.y-ix.y;
      r[2] = jx.z-ix.z;
      rnorm = gpu_dot3(r,r);
      rnorm = ucl_rsqrt(rnorm);
      rhat[0] = r[0]*rnorm;
      rhat[1] = r[1]*rnorm;
      rhat[2] = r[2]*rnorm;


      numtyp a2[9];       // Rotation matrix (lab->body)
      numtyp gamma2[9];   // A'*S^2*A
      numtyp4 jshape;

      jshape=shape[jtype];
      numtyp4 jshape2;
      jshape2.x=jshape.x*jshape.x;
      jshape2.y=jshape.y*jshape.y;
      jshape2.z=jshape.z*jshape.z;
      {
        numtyp aTs[9];    // A1'*S1^2
        gpu_quat_to_mat_trans(q,j,a2);
        gpu_transpose_times_diag3(a2,jshape2,aTs);
        gpu_times3(aTs,a2,gamma2);
      }

      numtyp temp[9], s[3], z1[3], z2[3], v1[3], v2[3];
      numtyp sigma12, sigma1, sigma2;
      gpu_plus3(gamma1,gamma2,temp);
      gpu_mldivide3(temp,rhat,s,err_flag);
      sigma12 = ucl_rsqrt((numtyp)0.5*gpu_dot3(s,rhat));
      gpu_times_column3(a1,rhat,z1);
      gpu_times_column3(a2,rhat,z2);
      v1[0] = z1[0]*ishape2.x;
      v1[1] = z1[1]*ishape2.y;
      v1[2] = z1[2]*ishape2.z;
      v2[0] = z2[0]/jshape2.x;
      v2[1] = z2[1]/jshape2.y;
      v2[2] = z2[2]/jshape2.z;
      sigma1 = ucl_sqrt(gpu_dot3(z1,v1));
      sigma2 = ucl_sqrt(gpu_dot3(z2,v2));

      numtyp H12[9];
      numtyp dH;
      H12[0] = gamma1[0]*sigma1+gamma2[0]*sigma2;
      H12[1] = gamma1[1]*sigma1+gamma2[1]*sigma2;
      H12[2] = gamma1[2]*sigma1+gamma2[2]*sigma2;
      H12[3] = gamma1[3]*sigma1+gamma2[3]*sigma2;
      H12[4] = gamma1[4]*sigma1+gamma2[4]*sigma2;
      H12[5] = gamma1[5]*sigma1+gamma2[5]*sigma2;
      H12[6] = gamma1[6]*sigma1+gamma2[6]*sigma2;
      H12[7] = gamma1[7]*sigma1+gamma2[7]*sigma2;
      H12[8] = gamma1[8]*sigma1+gamma2[8]*sigma2;
      dH=gpu_det3(H12);

      numtyp sigma1p2, sigma2p2, lambda, nu;
      sigma1p2 = sigma1*sigma1;
      sigma2p2 = sigma2*sigma2;
      numtyp jlshape = jshape.x*jshape.y*jshape.z;
      lambda = ilshape*sigma1p2 + jlshape*sigma2p2;


      sigma1=ucl_recip(sigma1);
      sigma2=ucl_recip(sigma2);

      nu = ucl_rsqrt((sigma1+sigma2)/dH);
      gpu_times3(aTe1,a1,temp);

      numtyp sigma, epsilon;
      int mtype=fast_mul(ntypes,itype)+jtype;
      sigma = sig_eps[mtype].x;
      epsilon = sig_eps[mtype].y*factor_lj;

      numtyp w[3], temp2[9];
      numtyp h12,eta,chi,sprod,sigh,tprod;
      numtyp aTe2[9];     // A'*E
      gpu_transpose_times_diag3(a2,well[jtype],aTe2);
      gpu_times3(aTe2,a2,temp2);
      gpu_plus3(temp,temp2,temp);
      gpu_mldivide3(temp,rhat,w,err_flag);
      h12 = ucl_recip(rnorm)-sigma12;
      eta = lambda/nu;
      chi = (numtyp)2.0*gpu_dot3(rhat,w);
      sprod = ilshape * jlshape;
      sigh = sigma/h12;
      tprod = eta*chi*sigh;

      numtyp stemp, Ua;
      stemp = h12*(numtyp)0.5;
      Ua = (ishape.x+stemp)*(ishape.y+stemp)*
           (ishape.z+stemp)*(jshape.x+stemp)*
           (jshape.y+stemp)*(jshape.z+stemp);
      Ua = ((numtyp)1.0+(numtyp)3.0*tprod)*sprod/Ua;
      Ua = epsilon*Ua/(numtyp)-36.0;

      numtyp Ur;
      stemp = h12/cr60;
      Ur = (ishape.x+stemp)*(ishape.y+stemp)*
           (ishape.z+stemp)*(jshape.x+stemp)*
           (jshape.y+stemp)*(jshape.z+stemp);
      Ur = ((numtyp)1.0+b_alpha*tprod)*sprod/Ur;
      numtyp sigh6=sigh*sigh*sigh;
      sigh6*=sigh6;
      Ur = epsilon*Ur*sigh6/(numtyp)2025.0;

      energy+=Ua+Ur;

      // force

      numtyp vsigma1[3], vsigma2[3], gsigma1[9], gsigma2[9];
      numtyp sec, sigma12p3, sigma1p3, sigma2p3;
      sec = sigma*eta*chi;
      sigma12p3 = sigma12*sigma12*sigma12;
      sigma1p3 = sigma1/sigma1p2;
      sigma2p3 = sigma2/sigma2p2;
      vsigma1[0] = -sigma1p3*v1[0];
      vsigma1[1] = -sigma1p3*v1[1];
      vsigma1[2] = -sigma1p3*v1[2];
      vsigma2[0] = -sigma2p3*v2[0];
      vsigma2[1] = -sigma2p3*v2[1];
      vsigma2[2] = -sigma2p3*v2[2];
      gsigma1[0] = -gamma1[0]*sigma1p2;
      gsigma1[1] = -gamma1[1]*sigma1p2;
      gsigma1[2] = -gamma1[2]*sigma1p2;
      gsigma1[3] = -gamma1[3]*sigma1p2;
      gsigma1[4] = -gamma1[4]*sigma1p2;
      gsigma1[5] = -gamma1[5]*sigma1p2;
      gsigma1[6] = -gamma1[6]*sigma1p2;
      gsigma1[7] = -gamma1[7]*sigma1p2;
      gsigma1[8] = -gamma1[8]*sigma1p2;
      gsigma2[0] = -gamma2[0]*sigma2p2;
      gsigma2[1] = -gamma2[1]*sigma2p2;
      gsigma2[2] = -gamma2[2]*sigma2p2;
      gsigma2[3] = -gamma2[3]*sigma2p2;
      gsigma2[4] = -gamma2[4]*sigma2p2;
      gsigma2[5] = -gamma2[5]*sigma2p2;
      gsigma2[6] = -gamma2[6]*sigma2p2;
      gsigma2[7] = -gamma2[7]*sigma2p2;
      gsigma2[8] = -gamma2[8]*sigma2p2;

      numtyp tsig1sig2, tdH, teta1, teta2;
      numtyp fourw[3], spr[3];
      tsig1sig2 = eta/((numtyp)2.0*(sigma1+sigma2));
      tdH = eta/((numtyp)2.0*dH);
      teta1 = (numtyp)2.0*eta/lambda;
      teta2 = teta1*jlshape/sigma2p3;
      teta1 = teta1*ilshape/sigma1p3;
      fourw[0] = (numtyp)4.0*w[0];
      fourw[1] = (numtyp)4.0*w[1];
      fourw[2] = (numtyp)4.0*w[2];
      spr[0] = (numtyp)0.5*sigma12p3*s[0];
      spr[1] = (numtyp)0.5*sigma12p3*s[1];
      spr[2] = (numtyp)0.5*sigma12p3*s[2];

      numtyp hsec, dspu, pbsu;
      stemp = ucl_recip(ishape.x*(numtyp)2.0+h12)+
              ucl_recip(ishape.y*(numtyp)2.0+h12)+
              ucl_recip(ishape.z*(numtyp)2.0+h12)+
              ucl_recip(jshape.x*(numtyp)2.0+h12)+
              ucl_recip(jshape.y*(numtyp)2.0+h12)+
              ucl_recip(jshape.z*(numtyp)2.0+h12);
      hsec = ucl_recip(h12+(numtyp)3.0*sec);
      dspu = ucl_recip(h12)-hsec+stemp;
      pbsu = (numtyp)3.0*sigma*hsec;

      numtyp dspr, pbsr;
      stemp = ucl_recip(ishape.x*cr60+h12)+
              ucl_recip(ishape.y*cr60+h12)+
              ucl_recip(ishape.z*cr60+h12)+
              ucl_recip(jshape.x*cr60+h12)+
              ucl_recip(jshape.y*cr60+h12)+
              ucl_recip(jshape.z*cr60+h12);
      hsec = ucl_recip(h12+b_alpha*sec);
      dspr = (numtyp)7.0/h12-hsec+stemp;
      pbsr = b_alpha*sigma*hsec;

      numtyp dH12[9];
      numtyp dUa, dUr, deta, dchi, ddH, dh12;
      numtyp dsigma1, dsigma2;

      #pragma unroll
      for (int i=0; i<3; i++) {
        numtyp u[3], u1[3], u2[3];
        u[0] = -rhat[i]*rhat[0];
        u[1] = -rhat[i]*rhat[1];
        u[2] = -rhat[i]*rhat[2];
        u[i] += (numtyp)1.0;
        u[0] *= rnorm;
        u[1] *= rnorm;
        u[2] *= rnorm;
        gpu_times_column3(a1,u,u1);
        gpu_times_column3(a2,u,u2);
        dsigma1=gpu_dot3(u1,vsigma1);
        dsigma2=gpu_dot3(u2,vsigma2);
        dH12[0] = dsigma1*gsigma1[0]+dsigma2*gsigma2[0];
        dH12[1] = dsigma1*gsigma1[1]+dsigma2*gsigma2[1];
        dH12[2] = dsigma1*gsigma1[2]+dsigma2*gsigma2[2];
        dH12[3] = dsigma1*gsigma1[3]+dsigma2*gsigma2[3];
        dH12[4] = dsigma1*gsigma1[4]+dsigma2*gsigma2[4];
        dH12[5] = dsigma1*gsigma1[5]+dsigma2*gsigma2[5];
        dH12[6] = dsigma1*gsigma1[6]+dsigma2*gsigma2[6];
        dH12[7] = dsigma1*gsigma1[7]+dsigma2*gsigma2[7];
        dH12[8] = dsigma1*gsigma1[8]+dsigma2*gsigma2[8];
        ddH = det_prime(H12,dH12);
        deta = (dsigma1+dsigma2)*tsig1sig2;
        deta -= ddH*tdH;
        deta -= dsigma1*teta1+dsigma2*teta2;
        dchi = gpu_dot3(u,fourw);
        dh12 = rhat[i]+gpu_dot3(u,spr);
        dUa = pbsu*(eta*dchi+deta*chi)-dh12*dspu;
        dUr = pbsr*(eta*dchi+deta*chi)-dh12*dspr;
        numtyp force=dUr*Ur+dUa*Ua;
        if (i==0) {
          f.x+=force;
          if (EVFLAG && vflag)
            virial[0]+=-r[0]*force;
        } else if (i==1) {
          f.y+=force;
          if (EVFLAG && vflag) {
            virial[1]+=-r[1]*force;
            virial[3]+=-r[0]*force;
          }
        } else {
          f.z+=force;
          if (EVFLAG && vflag) {
            virial[2]+=-r[2]*force;
            virial[4]+=-r[0]*force;
            virial[5]+=-r[1]*force;
          }
        }
      }

      // torque on i
      sigma1=ucl_recip(sigma1);

      numtyp fwae[3], p[3];
      gpu_row_times3(fourw,aTe1,fwae);

      {
        gpu_times_column3(lA1_0,rhat,p);
        dsigma1 = gpu_dot3(p,vsigma1);
        dH12[0] = lAsa1_0[0]*sigma1+dsigma1*gsigma1[0];
        dH12[1] = lAsa1_0[1]*sigma1+dsigma1*gsigma1[1];
        dH12[2] = lAsa1_0[2]*sigma1+dsigma1*gsigma1[2];
        dH12[3] = lAsa1_0[3]*sigma1+dsigma1*gsigma1[3];
        dH12[4] = lAsa1_0[4]*sigma1+dsigma1*gsigma1[4];
        dH12[5] = lAsa1_0[5]*sigma1+dsigma1*gsigma1[5];
        dH12[6] = lAsa1_0[6]*sigma1+dsigma1*gsigma1[6];
        dH12[7] = lAsa1_0[7]*sigma1+dsigma1*gsigma1[7];
        dH12[8] = lAsa1_0[8]*sigma1+dsigma1*gsigma1[8];
        ddH = det_prime(H12,dH12);
        deta = tsig1sig2*dsigma1-tdH*ddH;
        deta -= teta1*dsigma1;
        numtyp tempv[3];
        gpu_times_column3(lA1_0,w,tempv);
        dchi = -gpu_dot3(fwae,tempv);
        gpu_times_column3(lAtwo1_0,spr,tempv);
        dh12 = -gpu_dot3(s,tempv);

        dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
        dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
        tor.x -= (dUa*Ua+dUr*Ur);
      }

      {
        gpu_times_column3(lA1_1,rhat,p);
        dsigma1 = gpu_dot3(p,vsigma1);
        dH12[0] = lAsa1_1[0]*sigma1+dsigma1*gsigma1[0];
        dH12[1] = lAsa1_1[1]*sigma1+dsigma1*gsigma1[1];
        dH12[2] = lAsa1_1[2]*sigma1+dsigma1*gsigma1[2];
        dH12[3] = lAsa1_1[3]*sigma1+dsigma1*gsigma1[3];
        dH12[4] = lAsa1_1[4]*sigma1+dsigma1*gsigma1[4];
        dH12[5] = lAsa1_1[5]*sigma1+dsigma1*gsigma1[5];
        dH12[6] = lAsa1_1[6]*sigma1+dsigma1*gsigma1[6];
        dH12[7] = lAsa1_1[7]*sigma1+dsigma1*gsigma1[7];
        dH12[8] = lAsa1_1[8]*sigma1+dsigma1*gsigma1[8];
        ddH = det_prime(H12,dH12);
        deta = tsig1sig2*dsigma1-tdH*ddH;
        deta -= teta1*dsigma1;
        numtyp tempv[3];
        gpu_times_column3(lA1_1,w,tempv);
        dchi = -gpu_dot3(fwae,tempv);
        gpu_times_column3(lAtwo1_1,spr,tempv);
        dh12 = -gpu_dot3(s,tempv);

        dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
        dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
        tor.y -= (dUa*Ua+dUr*Ur);
      }

      {
        gpu_times_column3(lA1_2,rhat,p);
        dsigma1 = gpu_dot3(p,vsigma1);
        dH12[0] = lAsa1_2[0]*sigma1+dsigma1*gsigma1[0];
        dH12[1] = lAsa1_2[1]*sigma1+dsigma1*gsigma1[1];
        dH12[2] = lAsa1_2[2]*sigma1+dsigma1*gsigma1[2];
        dH12[3] = lAsa1_2[3]*sigma1+dsigma1*gsigma1[3];
        dH12[4] = lAsa1_2[4]*sigma1+dsigma1*gsigma1[4];
        dH12[5] = lAsa1_2[5]*sigma1+dsigma1*gsigma1[5];
        dH12[6] = lAsa1_2[6]*sigma1+dsigma1*gsigma1[6];
        dH12[7] = lAsa1_2[7]*sigma1+dsigma1*gsigma1[7];
        dH12[8] = lAsa1_2[8]*sigma1+dsigma1*gsigma1[8];
        ddH = det_prime(H12,dH12);
        deta = tsig1sig2*dsigma1-tdH*ddH;
        deta -= teta1*dsigma1;
        numtyp tempv[3];
        gpu_times_column3(lA1_2,w,tempv);
        dchi = -gpu_dot3(fwae,tempv);
        gpu_times_column3(lAtwo1_2,spr,tempv);
        dh12 = -gpu_dot3(s,tempv);

        dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
        dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
        tor.z -= (dUa*Ua+dUr*Ur);
      }

    } // for nbor
  } // if ii
  store_answers_t(f,tor,energy,virial,ii,astride,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv,inum);
}
