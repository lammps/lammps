// **************************************************************************
//                                gayberne_lj.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for Gay-Berne - Lennard-Jones potential acceleration
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_ellipsoid_extra.h"
#endif

__kernel void k_gayberne_sphere_ellipsoid(const __global numtyp4 *restrict x_,
                                          const __global numtyp4 *restrict q,
                                          const __global numtyp4 *restrict shape,
                                          const __global numtyp4 *restrict well,
                                          const __global numtyp *restrict gum,
                                          const __global numtyp2 *restrict sig_eps,
                                          const int ntypes,
                                          const __global numtyp *restrict lshape,
                                          const __global int *dev_nbor,
                                          const int stride,
                                          __global acctyp4 *restrict ans,
                                          __global acctyp *restrict engv,
                                          __global int *restrict err_flag,
                                          const int eflag, const int vflag,
                                          const int start, const int inum,
                                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  ii+=start;

  __local numtyp sp_lj[4];
  sp_lj[0]=gum[3];
  sp_lj[1]=gum[4];
  sp_lj[2]=gum[5];
  sp_lj[3]=gum[6];

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info_e(dev_nbor,stride,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex);
    int itype=ix.w;

    numtyp oner=shape[itype].x;
    numtyp one_well=well[itype].x;

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_nbor[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex);
      int jtype=jx.w;

      // Compute r12
      numtyp r12[3];
      r12[0] = jx.x-ix.x;
      r12[1] = jx.y-ix.y;
      r12[2] = jx.z-ix.z;
      numtyp ir = gpu_dot3(r12,r12);

      ir = ucl_rsqrt(ir);
      numtyp r = ucl_recip(ir);

      numtyp r12hat[3];
      r12hat[0]=r12[0]*ir;
      r12hat[1]=r12[1]*ir;
      r12hat[2]=r12[2]*ir;

      numtyp a2[9];
      gpu_quat_to_mat_trans(q,j,a2);

      numtyp u_r, dUr[3], eta;
      { // Compute U_r, dUr, eta, and teta
        // Compute g12
        numtyp g12[9];
        {
          {
            numtyp g2[9];
            gpu_diag_times3(shape[jtype],a2,g12);
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
              sigma12 = ucl_rsqrt((numtyp)0.5*sigma12);
              h12 = r-sigma12;

              // -- kappa is now ok
              kappa[0]*=r;
              kappa[1]*=r;
              kappa[2]*=r;

              int mtype=fast_mul(ntypes,itype)+jtype;
              numtyp sigma = sig_eps[mtype].x;
              numtyp epsilon = sig_eps[mtype].y;
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
          eta = (numtyp)2.0*lshape[itype]*lshape[jtype];
          numtyp det_g12 = gpu_det3(g12);
          eta = ucl_powr(eta/det_g12,gum[1]);
        }
      }

      numtyp chi, dchi[3];
      { // Compute chi and dchi

        // Compute b12
        numtyp b12[9];
        {
          numtyp b2[9];
          gpu_diag_times3(well[jtype],a2,b12);
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
        chi = ucl_powr(chi*(numtyp)2.0,gum[2]);

        // -- iota is now ok
        iota[0]*=r;
        iota[1]*=r;
        iota[2]*=r;

        numtyp temp1 = gpu_dot3(iota,r12hat);
        numtyp temp2 = (numtyp)-4.0*ir*ir*gum[2]*ucl_powr(chi,(gum[2]-(numtyp)1.0)/
                                                     gum[2]);
        dchi[0] = temp2*(iota[0]-temp1*r12hat[0]);
        dchi[1] = temp2*(iota[1]-temp1*r12hat[1]);
        dchi[2] = temp2*(iota[2]-temp1*r12hat[2]);
      }

      numtyp temp2 = factor_lj*eta*chi;
      if (eflag>0)
        energy+=u_r*temp2;
      numtyp temp1 = -eta*u_r*factor_lj;
      if (vflag>0) {
        r12[0]*=-1;
        r12[1]*=-1;
        r12[2]*=-1;
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
    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_gayberne_lj(const __global numtyp4 *restrict x_,
                            const __global numtyp4 *restrict lj1,
                            const __global numtyp4 *restrict lj3,
                            const int lj_types,
                            const __global numtyp *restrict gum,
                            const int stride,
                            const __global int *dev_ij,
                            __global acctyp4 *restrict ans,
                            __global acctyp *restrict engv,
                            __global int *restrict err_flag,
                            const int eflag, const int vflag, const int start,
                            const int inum, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  ii+=start;

  __local numtyp sp_lj[4];
  sp_lj[0]=gum[3];
  sp_lj[1]=gum[4];
  sp_lj[2]=gum[5];
  sp_lj[3]=gum[6];

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info_e(dev_ij,stride,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex);
    int itype=ix.w;

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_ij[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex);
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp r2inv = delx*delx+dely*dely+delz*delz;

      int ii=itype*lj_types+jtype;
      if (r2inv<lj1[ii].z && lj1[ii].w==SPHERE_SPHERE) {
        r2inv=ucl_recip(r2inv);
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = r2inv*r6inv*(lj1[ii].x*r6inv-lj1[ii].y);
        force*=factor_lj;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=r6inv*(lj3[ii].x*r6inv-lj3[ii].y);
          energy+=factor_lj*(e-lj3[ii].z);
        }
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
    acc_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
  } // if ii
}

__kernel void k_gayberne_lj_fast(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict lj1_in,
                                 const __global numtyp4 *restrict lj3_in,
                                 const __global numtyp *restrict gum,
                                 const int stride,
                                 const __global int *dev_ij,
                                 __global acctyp4 *restrict ans,
                                 __global acctyp *restrict engv,
                                 __global int *restrict err_flag,
                                 const int eflag, const int vflag,
                                 const int start, const int inum,
                                 const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  ii+=start;

  __local numtyp sp_lj[4];
  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (tid<4)
    sp_lj[tid]=gum[tid+3];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (eflag>0)
      lj3[tid]=lj3_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info_e(dev_ij,stride,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex);
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_ij[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex);
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp r2inv = delx*delx+dely*dely+delz*delz;

      if (r2inv<lj1[mtype].z && lj1[mtype].w==SPHERE_SPHERE) {
        r2inv=ucl_recip(r2inv);
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = factor_lj*r2inv*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          energy+=factor_lj*(e-lj3[mtype].z);
        }
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
    acc_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
  } // if ii
}

