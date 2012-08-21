// **************************************************************************
//                                coul_dsf.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the coul/dsf pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 8/15/2012
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL

#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float> q_tex;
#else
texture<int4,1> pos_tex;
texture<int2> q_tex;
#endif

#else
#define pos_tex x_
#define q_tex q_
#endif

#define MY_PIS (acctyp)1.77245385090551602729

__kernel void k_coul_dsf(__global numtyp4 *x_, const int lj_types, 
                         __global numtyp *sp_lj_in, __global int *dev_nbor, 
                         __global int *dev_packed, __global acctyp4 *ans,
                         __global acctyp *engv, const int eflag,
                         const int vflag, const int inum,
                         const int nbor_pitch, __global numtyp *q_ ,
                         const numtyp cut_coulsq, const numtyp qqrd2e,
                         const numtyp e_shift, const numtyp f_shift, 
                         const numtyp alpha, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;

      numtyp factor_coul, r, prefactor, erfcc;
      factor_coul = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq < cut_coulsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force;

        r = ucl_sqrt(rsq);
        fetch(prefactor,j,q_tex);
        prefactor *= factor_coul * qqrd2e*qtmp/r;
        numtyp erfcd = ucl_exp(-alpha*alpha*rsq);
        numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*alpha*r);
        erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
        forcecoul = prefactor * (erfcc + 2.0*alpha/MY_PIS*r*erfcd + 
                                 rsq*f_shift);
        
        force = forcecoul * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          if (rsq < cut_coulsq) {
            numtyp e=prefactor*(erfcc-r*e_shift-rsq*f_shift);
            e_coul += e;
          }
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

__kernel void k_coul_dsf_fast(__global numtyp4 *x_, __global numtyp* sp_lj_in,
                              __global int *dev_nbor, __global int *dev_packed,
                              __global acctyp4 *ans, __global acctyp *engv, 
                              const int eflag, const int vflag, const int inum, 
                              const int nbor_pitch, __global numtyp *q_,
                              const numtyp cut_coulsq, const numtyp qqrd2e,
                              const numtyp e_shift, const numtyp f_shift, 
                              const numtyp alpha, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  
  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __syncthreads();
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
 
    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;

      numtyp factor_coul, r, prefactor, erfcc;
      factor_coul = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq < cut_coulsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force;

        r = ucl_sqrt(rsq);
        fetch(prefactor,j,q_tex);
        prefactor *= factor_coul * qqrd2e*qtmp/r;
        numtyp erfcd = ucl_exp(-alpha*alpha*rsq);
        numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*alpha*r);
        erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
        forcecoul = prefactor * (erfcc + 2.0*alpha/MY_PIS*r*erfcd + 
                                   rsq*f_shift);
        
        force = forcecoul * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          if (rsq < cut_coulsq) {
            numtyp e=prefactor*(erfcc-r*e_shift-rsq*f_shift);
            e_coul += e;
          }
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

