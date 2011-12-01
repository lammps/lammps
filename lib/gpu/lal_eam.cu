// **************************************************************************
//                             lal_eam.cu
//                             -------------------
//                     Trung Dac Nguyen, W. Michael Brown (ORNL)
//
//  Device code for acceleration of the eam pair style
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
texture<float> fp_tex;
#ifndef _DOUBLE_DOUBLE
ucl_inline float4 fetch_pos(const int& i, const float4 *pos)
  { return tex1Dfetch(pos_tex, i); }
ucl_inline float fetch_fp(const int& i, const float *fp) 
  { return tex1Dfetch(fp_tex, i); }
#endif
#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define store_answers_eam(f, energy, e_coul, virial, ii, inum, tid,         \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    red_acc[4][tid]=e_coul;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<5; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
    e_coul=red_acc[4][tid];                                                 \
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
    engv+=ii;                                                               \
    if (eflag>0) {                                                          \
      *engv+=energy;                                                        \
      engv+=inum;                                                           \
      *engv+=e_coul;                                                        \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *engv=virial[i];                                                    \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

__kernel void kernel_energy(__global numtyp4 *x_, 
                    __global numtyp2 *type2rhor_z2r, __global numtyp *type2frho,
                    __global numtyp *rhor_spline, __global numtyp *frho_spline,
                    __global int *dev_nbor, __global int *dev_packed,
                    __global acctyp *fp_, 
                    __global acctyp *engv, const int eflag, 
                    const int vflag, const int inum, 
                    const int nbor_pitch,
                    const int ntypes, const numtyp cutforcesq, 
                    const numtyp rdr, const numtyp rdrho,
                    const int nrho, const int nr,
                    const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp rho = (acctyp)0;
  acctyp energy = (acctyp)0;
   
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int itype=ix.w;
    
    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
      
      if (rsq<cutforcesq) {
        numtyp p = ucl_sqrt(rsq)*rdr + (numtyp)1.0;
        int m=p;
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,(numtyp)1.0);
        
        int mtype = jtype*ntypes+itype;
        int index = type2rhor_z2r[mtype].x*(nr+1)*7+m*7;
        numtyp coeff3 = rhor_spline[index+3];
        numtyp coeff4 = rhor_spline[index+4];
        numtyp coeff5 = rhor_spline[index+5];
        numtyp coeff6 = rhor_spline[index+6];
        rho += ((coeff3*p + coeff4)*p + coeff5)*p + coeff6;
      }
    } // for nbor
    
    // reduce to get the density at atom ii
    
    if (t_per_atom>1) {
      __local acctyp red_acc[BLOCK_PAIR];                               
      red_acc[tid]=rho;                                                 
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                      
        if (offset < s)                                        
          red_acc[tid] += red_acc[tid+s];
      }                                                                    
      rho=red_acc[tid];                                                 
    }
    
    // calculate the embedded force for ii
    if (offset==0) {
      numtyp p = rho*rdrho + (numtyp)1.0;
      int m=p;
      m = MAX(1,MIN(m,nrho-1));
      p -= m;
      p = MIN(p,(numtyp)1.0);
      
      int index = type2frho[itype]*(nr+1)*7+m*7;
      numtyp coeff0 = frho_spline[index+0];
      numtyp coeff1 = frho_spline[index+1];
      numtyp coeff2 = frho_spline[index+2];
      numtyp fp = (coeff0*p + coeff1)*p + coeff2;
      fp_[ii]=fp;                       
      
      engv+=ii;  
      if (eflag>0) {
        numtyp coeff3 = frho_spline[index+3];
        numtyp coeff4 = frho_spline[index+4];
        numtyp coeff5 = frho_spline[index+5];
        numtyp coeff6 = frho_spline[index+6];
        energy = ((coeff3*p + coeff4)*p + coeff5)*p + coeff6;
        *engv=(acctyp)2.0*energy;
      }
    }
  } // if ii
}


__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp *fp_,
                          __global numtyp2 *type2rhor_z2r,
                          __global numtyp *rhor_spline, __global numtyp *z2r_spline,
                          __global int *dev_nbor, 
                          __global int *dev_packed, __global acctyp4 *ans,
                          __global acctyp *engv, const int eflag, 
                          const int vflag, const int inum, const int nbor_pitch,
                          const int ntypes, const numtyp cutforcesq, 
                          const numtyp rdr, const int nr,
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
    numtyp ifp=fetch_fp(i,fp_);  //fp_[i];
    int itype=ix.w;

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      numtyp jfp=fetch_fp(j,fp_); //fp_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
      
      if (rsq<cutforcesq) {
        numtyp r = ucl_sqrt(rsq);
        numtyp p = r*rdr + (numtyp)1.0;
        int m=p;
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
    store_answers_eam(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
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
                          const numtyp rdr, const int nr,
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
    numtyp ifp=fetch_fp(i,fp_); //fp_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      numtyp jfp=fetch_fp(j,fp_); //fp_[j];
      int jtype=jx.w;
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      if (rsq<cutforcesq) {
        numtyp r = ucl_sqrt(rsq);
        numtyp p = r*rdr + (numtyp)1.0;
        int m=p;
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
    store_answers_eam(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

