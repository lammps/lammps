// **************************************************************************
//                             lal_eam.cu
//                             -------------------
//                     W. Michael Brown, Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the lj/cut/coul/fsww/cut pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
texture<float4> pos_tex;
texture<float> q_tex;
#ifndef _DOUBLE_DOUBLE
ucl_inline float4 fetch_pos(const int& i, const float4 *pos)
  { return tex1Dfetch(pos_tex, i); }
ucl_inline float fetch_q(const int& i, const float *q) 
  { return tex1Dfetch(q_tex, i); }
#endif
#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))

__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp *fp_,
                          __global numtyp2 *type2rhor_z2r,
                          __global numtyp *rhor_spline, __global numtyp *z2r_spline,
                          __global int *dev_nbor, 
                          __global int *dev_packed, __global acctyp4 *ans,
                          __global acctyp *engv, const int eflag, 
                          const int vflag, const int inum, const int nbor_pitch,
                          const int ntypes, const numtyp cutforcesq, 
                          const numtyp rdr, const int nrhor, const int nz2r, const int nr,
                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    numtyp ifp=fetch_q(i,fp_);  //fp_[i];
    int itype=ix.w;

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      numtyp jfp=fetch_q(j,fp_); //fp_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
      
      if (rsq<cutforcesq) {
        numtyp r = ucl_sqrt(rsq);
        numtyp p = r*rdr + (numtyp)1.0;
        int m=__float2int_rn(p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,(numtyp)1.0);
        
        int mtype,index;
        numtyp coeff0,coeff1,coeff2,coeff3,coeff4,coeff5,coeff6;
        
        mtype = itype*ntypes+jtype;
        index = type2rhor_z2r[mtype].x*(nr+1)*7+m*7;
        coeff0 = rhor_spline[index+0];
        coeff1 = rhor_spline[index+1];
        coeff2 = rhor_spline[index+2];
        numtyp rhoip = (coeff0*p + coeff1)*p + coeff2;
        
        mtype = jtype*ntypes+itype;
        index = type2rhor_z2r[mtype].x*(nr+1)*7+m*7;
        coeff0 = rhor_spline[index+0];
        coeff1 = rhor_spline[index+1];
        coeff2 = rhor_spline[index+2];
        numtyp rhojp = (coeff0*p + coeff1)*p + coeff2;
        
        mtype = itype*ntypes+jtype;
        index = type2rhor_z2r[mtype].y*(nr+1)*7+m*7;
        coeff0 = z2r_spline[index+0];
        coeff1 = z2r_spline[index+1];
        coeff2 = z2r_spline[index+2];
        coeff3 = z2r_spline[index+3];
        coeff4 = z2r_spline[index+4];
        coeff5 = z2r_spline[index+5];
        coeff6 = z2r_spline[index+6];
        
        numtyp z2p = (coeff0*p + coeff1)*p + coeff2;
        numtyp z2 = ((coeff3*p + coeff4)*p + coeff5)*p + coeff6;

        numtyp recip = (numtyp)1.0/r;
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip = ifp*rhojp + jfp*rhoip + phip;
        numtyp force = -psip*recip;
        
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          energy += phi;
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii

}

__kernel void kernel_pair_fast(__global numtyp4 *x_, __global numtyp *fp_,
                          __global numtyp2 *type2rhor_z2r,
                          __global numtyp *rhor_spline, __global numtyp *z2r_spline,
                          __global int *dev_nbor, __global int *dev_packed, 
                          __global acctyp4 *ans, __global acctyp *engv, 
                          const int eflag, const int vflag, const int inum, 
                          const int nbor_pitch,
                          const numtyp cutforcesq, 
                          const numtyp rdr, const int nrhor, const int nz2r, const int nr,
                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
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

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    numtyp ifp=fetch_q(i,fp_); //fp_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      numtyp jfp=fetch_q(j,fp_); //fp_[j];
      int jtype=jx.w;
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      if (rsq<cutforcesq) {
        numtyp r = ucl_sqrt(rsq);
        numtyp p = r*rdr + (numtyp)1.0;
        int m=__float2int_rn(p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,(numtyp)1.0);
        
        numtyp coeff0,coeff1,coeff2,coeff3,coeff4,coeff5,coeff6;
        int mtype,index;
        
        mtype = itype+jx.w;
        index = type2rhor_z2r[mtype].x*(nr+1)*7+m*7;
        coeff0 = rhor_spline[index+0];
        coeff1 = rhor_spline[index+1];
        coeff2 = rhor_spline[index+2];
        numtyp rhoip = (coeff0*p + coeff1)*p + coeff2;
        
        mtype = jtype+ix.w;
        index = type2rhor_z2r[mtype].x*(nr+1)*7+m*7;
        coeff0 = rhor_spline[index+0];
        coeff1 = rhor_spline[index+1];
        coeff2 = rhor_spline[index+2];
        numtyp rhojp = (coeff0*p + coeff1)*p + coeff2;
        
        mtype = itype+jx.w;
        index = type2rhor_z2r[mtype].y*(nr+1)*7+m*7;
        coeff0 = z2r_spline[index+0];
        coeff1 = z2r_spline[index+1];
        coeff2 = z2r_spline[index+2];
        coeff3 = z2r_spline[index+3];
        coeff4 = z2r_spline[index+4];
        coeff5 = z2r_spline[index+5];
        coeff6 = z2r_spline[index+6];
        
        numtyp z2p = (coeff0*p + coeff1)*p + coeff2;
        numtyp z2 = ((coeff3*p + coeff4)*p + coeff5)*p + coeff6;

        numtyp recip = 1.0/r;
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip = ifp*rhojp + jfp*rhoip + phip;
        numtyp force = -psip*recip;
        
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          energy += phi;
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

