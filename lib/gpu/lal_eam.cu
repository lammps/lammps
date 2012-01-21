// **************************************************************************
//                                   eam.cu
//                             -------------------
//                   Trung Dac Nguyen, W. Michael Brown (ORNL)
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

texture<float4> rhor_sp1_tex;
texture<float4> rhor_sp2_tex;
texture<float4> frho_sp1_tex;
texture<float4> frho_sp2_tex;
texture<float4> z2r_sp1_tex;
texture<float4> z2r_sp2_tex;

#ifdef _DOUBLE_DOUBLE
ucl_inline double4 fetch_rhor_sp1(const int& i, const double4 *rhor_spline1) { 
  return rhor_spline1[i]; 
}
ucl_inline double4 fetch_rhor_sp2(const int& i, const double4 *rhor_spline2) { 
  return rhor_spline2[i]; 
}
ucl_inline double4 fetch_frho_sp1(const int& i, const double4 *frho_spline1) { 
  return frho_spline1[i]; 
}
ucl_inline double4 fetch_frho_sp2(const int& i, const double4 *frho_spline2) { 
  return frho_spline2[i]; 
}
ucl_inline double4 fetch_z2r_sp1(const int& i, const double4 *z2r_spline1) { 
  return z2r_spline1[i]; 
}
ucl_inline double4 fetch_z2r_sp2(const int& i, const double4 *z2r_spline2) { 
  return z2r_spline2[i]; 
}
#endif

#ifndef _DOUBLE_DOUBLE
ucl_inline float4 fetch_pos(const int& i, const float4 *pos)
  { return tex1Dfetch(pos_tex, i); }
ucl_inline float fetch_q(const int& i, const float *fp) 
  { return tex1Dfetch(fp_tex, i); }

ucl_inline float4 fetch_rhor_sp1(const int& i, const float4 *rhor_spline1) 
  { return tex1Dfetch(rhor_sp1_tex, i); }
ucl_inline float4 fetch_rhor_sp2(const int& i, const float4 *rhor_spline2) 
  { return tex1Dfetch(rhor_sp2_tex, i); }
ucl_inline float4 fetch_frho_sp1(const int& i, const float4 *frho_spline1) 
  { return tex1Dfetch(frho_sp1_tex, i); }
ucl_inline float4 fetch_frho_sp2(const int& i, const float4 *frho_spline2) 
  { return tex1Dfetch(frho_sp2_tex, i); }
ucl_inline float4 fetch_z2r_sp1(const int& i, const float4 *z2r_spline1) 
  { return tex1Dfetch(z2r_sp1_tex, i); }
ucl_inline float4 fetch_z2r_sp2(const int& i, const float4 *z2r_spline2) 
  { return tex1Dfetch(z2r_sp2_tex, i); }
#endif

#else // OPENCL

#define fetch_q(i,y) fp_[i]
#define fetch_rhor_sp1(i,y) rhor_spline1[i]
#define fetch_rhor_sp2(i,y) rhor_spline2[i]
#define fetch_frho_sp1(i,y) frho_spline1[i]
#define fetch_frho_sp2(i,y) frho_spline2[i]
#define fetch_z2r_sp1(i,y) z2r_spline1[i] 
#define fetch_z2r_sp2(i,y) z2r_spline2[i]

#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#define store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,           \
                        eflag,vflag,engv,rdrho,nrho,i)                      \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[BLOCK_PAIR];                                     \
    red_acc[tid]=rho;                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s)                                                       \
         red_acc[tid] += red_acc[tid+s];                                    \
      }                                                                     \
      rho=red_acc[tid];                                                     \
  }                                                                         \
  if (offset==0) {                                                          \
    numtyp p = rho*rdrho + (numtyp)1.0;                                     \
    int m=p;                                                                \
    m = MAX(1,MIN(m,nrho-1));                                               \
    p -= m;                                                                 \
    p = MIN(p,(numtyp)1.0);                                                 \
    int index = type2frho[itype]*(nrho+1)+m;                                \
    numtyp4 coeff = fetch_frho_sp1(index, frho_spline1);                    \
    numtyp fp = (coeff.x*p + coeff.y)*p + coeff.z;                          \
    fp_[i]=fp;                                                              \
    if (eflag>0) {                                                          \
      coeff = fetch_frho_sp2(index, frho_spline2);                          \
      energy = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;             \
      engv[ii]=(acctyp)2.0*energy;                                          \
    }                                                                       \
  }

#define store_answers_eam(f, energy, virial, ii, inum, tid, t_per_atom,     \
                      offset, elag, vflag, ans, engv)                       \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<4; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
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
    if (eflag>0) {                                                          \
      engv[ii]+=energy;                                                     \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ii]=virial[i];                                                 \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

__kernel void kernel_energy(__global numtyp4 *x_, __global int2 *type2rhor_z2r,
                            __global int *type2frho, 
                            __global numtyp4 *rhor_spline2, 
                            __global numtyp4 *frho_spline1,
                            __global numtyp4 *frho_spline2,
                            __global int *dev_nbor, __global int *dev_packed,
                            __global numtyp *fp_, __global acctyp *engv, 
                            const int eflag, const int inum, 
                            const int nbor_pitch, const int ntypes, 
                            const numtyp cutforcesq, const numtyp rdr, 
                            const numtyp rdrho, const int nrho, const int nr,
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
        int index = type2rhor_z2r[mtype].x*(nr+1)+m;
        numtyp4 coeff = fetch_rhor_sp2(index, rhor_spline2);
        rho += ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      }
    } // for nbor
    
    store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,
        eflag,vflag,engv,rdrho,nrho,i);
  } // if ii
}

__kernel void kernel_energy_fast(__global numtyp4 *x_, 
                                 __global int2 *type2rhor_z2r_in,
                                 __global int *type2frho_in, 
                                 __global numtyp4 *rhor_spline2, 
                                 __global numtyp4 *frho_spline1,
                                 __global numtyp4 *frho_spline2,
                                 __global int *dev_nbor, 
                                 __global int *dev_packed, __global numtyp *fp_, 
                                 __global acctyp *engv, const int eflag, 
                                 const int inum, const int nbor_pitch,
                                 const int ntypes, const numtyp cutforcesq, 
                                 const numtyp rdr, const numtyp rdrho,
                                 const int nrho, const int nr, 
                                 const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local int2 type2rhor_z2r[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local int type2frho[MAX_SHARED_TYPES];

  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    type2rhor_z2r[tid]=type2rhor_z2r_in[tid];
  }
  
  if (tid<MAX_SHARED_TYPES) {
    type2frho[tid]=type2frho_in[tid];
  }

  acctyp rho = (acctyp)0;
  acctyp energy = (acctyp)0;
  
  __syncthreads(); 

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
        
        int jtype=fast_mul((int)MAX_SHARED_TYPES,jx.w);
        int mtype = jtype+itype;
        int index = type2rhor_z2r[mtype].x*(nr+1)+m;
        numtyp4 coeff = fetch_rhor_sp2(index, rhor_spline2);
        rho += ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      }
    } // for nbor
    
    store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,
                    eflag,vflag,engv,rdrho,nrho,i);
  } // if ii
}

__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp *fp_,
                          __global int2 *type2rhor_z2r,
                          __global numtyp4 *rhor_spline1, 
                          __global numtyp4 *z2r_spline1,
                          __global numtyp4 *z2r_spline2,
                          __global int *dev_nbor, __global int *dev_packed, 
                          __global acctyp4 *ans, __global acctyp *engv, 
                          const int eflag, const int vflag, 
                          const int inum, const int nbor_pitch,
                          const int ntypes, const numtyp cutforcesq, 
                          const numtyp rdr, const int nr,
                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  acctyp energy=(acctyp)0;
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
        numtyp4 coeff;

        mtype = itype*ntypes+jtype;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        coeff = fetch_rhor_sp1(index, rhor_spline1); 
        numtyp rhoip = (coeff.x*p + coeff.y)*p + coeff.z;

        mtype = jtype*ntypes+itype;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        coeff = fetch_rhor_sp1(index, rhor_spline1); 
        numtyp rhojp = (coeff.x*p + coeff.y)*p + coeff.z;
              
        mtype = itype*ntypes+jtype;
        index = type2rhor_z2r[mtype].y*(nr+1)+m;
        coeff = fetch_z2r_sp1(index, z2r_spline1);
        numtyp z2p = (coeff.x*p + coeff.y)*p + coeff.z;
        coeff = fetch_z2r_sp2(index, z2r_spline2);
        numtyp z2 = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
        
        numtyp recip = ucl_recip(r);
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip = ifp*rhojp + fetch_q(j,fp_)*rhoip + phip; 
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
    store_answers_eam(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii

}

__kernel void kernel_pair_fast(__global numtyp4 *x_, __global numtyp *fp_,
                          __global int2 *type2rhor_z2r_in,
                          __global numtyp4 *rhor_spline1, 
                          __global numtyp4 *z2r_spline1,
                          __global numtyp4 *z2r_spline2,
                          __global int *dev_nbor, __global int *dev_packed, 
                          __global acctyp4 *ans, __global acctyp *engv, 
                          const int eflag, const int vflag, const int inum, 
                          const int nbor_pitch,
                          const numtyp cutforcesq, 
                          const numtyp rdr, const int nr,
                          const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local int2 type2rhor_z2r[MAX_SHARED_TYPES*MAX_SHARED_TYPES];

  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    type2rhor_z2r[tid]=type2rhor_z2r_in[tid];
  }

  acctyp energy=(acctyp)0;
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

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    numtyp ifp=fetch_q(i,fp_); //fp_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jw=jx.w;
      int jtype=fast_mul((int)MAX_SHARED_TYPES,jw);
      
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
        
        numtyp4 coeff;
        int mtype,index;
        
        mtype = itype+jw;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        coeff = fetch_rhor_sp1(index, rhor_spline1); 
        numtyp rhoip = (coeff.x*p + coeff.y)*p + coeff.z;
        
        mtype = jtype+iw;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        coeff = fetch_rhor_sp1(index, rhor_spline1); 
        numtyp rhojp = (coeff.x*p + coeff.y)*p + coeff.z;
        
        mtype = itype+jw;
        index = type2rhor_z2r[mtype].y*(nr+1)+m;
        coeff = fetch_z2r_sp1(index, z2r_spline1);
        numtyp z2p = (coeff.x*p + coeff.y)*p + coeff.z;
        coeff = fetch_z2r_sp2(index, z2r_spline2);
        numtyp z2 = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      
        numtyp recip = ucl_recip(r);
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip = ifp*rhojp + fetch_q(j,fp_)*rhoip + phip;
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
    store_answers_eam(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

