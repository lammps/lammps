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

#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float> fp_tex;
texture<float4> rhor_sp1_tex;
texture<float4> rhor_sp2_tex;
texture<float4> frho_sp1_tex;
texture<float4> frho_sp2_tex;
texture<float4> z2r_sp1_tex;
texture<float4> z2r_sp2_tex;
#else
texture<int4> pos_tex;
texture<int2> fp_tex;
texture<int4> rhor_sp1_tex;
texture<int4> rhor_sp2_tex;
texture<int4> frho_sp1_tex;
texture<int4> frho_sp2_tex;
texture<int4> z2r_sp1_tex;
texture<int4> z2r_sp2_tex;
#endif

#else

#define pos_tex x_
#define fp_tex fp_
#define rhor_sp1_tex rhor_spline1
#define rhor_sp2_tex rhor_spline2
#define frho_sp1_tex frho_spline1
#define frho_sp2_tex frho_spline2
#define z2r_sp1_tex z2r_spline1
#define z2r_sp2_tex z2r_spline2

#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#if (ARCH < 300)

#define store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,           \
                        eflag,vflag,engv,rdrho,nrho,i,rhomax)               \
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
    numtyp4 coeff; fetch4(coeff,index,frho_sp1_tex);                        \
    numtyp fp = (coeff.x*p + coeff.y)*p + coeff.z;                          \
    fp_[i]=fp;                                                              \
    if (eflag>0) {                                                          \
      fetch4(coeff,index,frho_sp2_tex);                                     \
      energy = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;             \
      if (rho > rhomax) energy += fp*(rho-rhomax);                          \
      engv[ii]=energy;                                                      \
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
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#else

#define store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,           \
                        eflag,vflag,engv,rdrho,nrho,i,rhomax)               \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1)                           \
        rho += shfl_xor(rho, s, t_per_atom);                                \
  }                                                                         \
  if (offset==0) {                                                          \
    numtyp p = rho*rdrho + (numtyp)1.0;                                     \
    int m=p;                                                                \
    m = MAX(1,MIN(m,nrho-1));                                               \
    p -= m;                                                                 \
    p = MIN(p,(numtyp)1.0);                                                 \
    int index = type2frho[itype]*(nrho+1)+m;                                \
    numtyp4 coeff; fetch4(coeff,index,frho_sp1_tex);                        \
    numtyp fp = (coeff.x*p + coeff.y)*p + coeff.z;                          \
    fp_[i]=fp;                                                              \
    if (eflag>0) {                                                          \
      fetch4(coeff,index,frho_sp2_tex);                                     \
      energy = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;             \
      if (rho > rhomax) energy += fp*(rho-rhomax);                          \
      engv[ii]=energy;                                          \
    }                                                                       \
  }

#define store_answers_eam(f, energy, virial, ii, inum, tid, t_per_atom,     \
                          offset, eflag, vflag, ans, engv)                  \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        energy += shfl_xor(energy, s, t_per_atom);                          \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#endif

__kernel void k_energy(const __global numtyp4 *restrict x_, 
                       const __global int2 *restrict type2rhor_z2r,
                       const __global int *restrict type2frho, 
                       const __global numtyp4 *restrict rhor_spline2, 
                       const __global numtyp4 *restrict frho_spline1,
                       const __global numtyp4 *restrict frho_spline2,
                       const __global int *dev_nbor, 
                       const __global int *dev_packed,
                       __global numtyp *restrict fp_, 
                       __global acctyp *restrict engv, 
                       const int eflag, const int inum, const int nbor_pitch,
                       const int ntypes,  const numtyp cutforcesq, 
                       const numtyp rdr, const numtyp rdrho, 
                       const numtyp rhomax, const int nrho,
                       const int nr, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp rho = (acctyp)0;
  acctyp energy = (acctyp)0;
   
  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
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
        numtyp4 coeff; fetch4(coeff,index,rhor_sp2_tex);
        rho += ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      }
    } // for nbor
    
    store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,
        eflag,vflag,engv,rdrho,nrho,i,rhomax);
  } // if ii
}

__kernel void k_energy_fast(const __global numtyp4 *restrict x_, 
                            const __global int2 *restrict type2rhor_z2r_in,
                            const __global int *restrict type2frho_in, 
                            const __global numtyp4 *restrict rhor_spline2, 
                            const __global numtyp4 *restrict frho_spline1,
                            const __global numtyp4 *restrict frho_spline2,
                            const __global int *dev_nbor, 
                            const __global int *dev_packed, 
                            __global numtyp *restrict fp_, 
                            __global acctyp *restrict engv, 
                            const int eflag,  const int inum, 
                            const int nbor_pitch, const int ntypes, 
                            const numtyp cutforcesq,  const numtyp rdr,
                            const numtyp rdrho, const numtyp rhomax, 
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
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

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
        numtyp4 coeff; fetch4(coeff,index,rhor_sp2_tex);
        rho += ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      }
    } // for nbor
    
    store_energy_fp(rho,energy,ii,inum,tid,t_per_atom,offset,
                    eflag,vflag,engv,rdrho,nrho,i,rhomax);
  } // if ii
}

__kernel void k_eam(const __global numtyp4 *restrict x_, 
                    const __global numtyp *fp_,
                    const __global int2 *type2rhor_z2r,
                    const __global numtyp4 *rhor_spline1, 
                    const __global numtyp4 *z2r_spline1,
                    const __global numtyp4 *z2r_spline2,
                    const __global int *dev_nbor,
                    const __global int *dev_packed, 
                    __global acctyp4 *ans,
                    __global acctyp *engv, 
                    const int eflag, const int vflag,  const int inum,
                    const int nbor_pitch, const int ntypes, 
                    const numtyp cutforcesq,  const numtyp rdr, const int nr,
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
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);
  
    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp ifp; fetch(ifp,i,fp_tex);  //fp_[i];
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
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
        fetch4(coeff,index,rhor_sp1_tex);
        numtyp rhoip = (coeff.x*p + coeff.y)*p + coeff.z;

        mtype = jtype*ntypes+itype;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        fetch4(coeff,index,rhor_sp1_tex);
        numtyp rhojp = (coeff.x*p + coeff.y)*p + coeff.z;
              
        mtype = itype*ntypes+jtype;
        index = type2rhor_z2r[mtype].y*(nr+1)+m;
        fetch4(coeff,index,z2r_sp1_tex);
        numtyp z2p = (coeff.x*p + coeff.y)*p + coeff.z;
        fetch4(coeff,index,z2r_sp2_tex);
        numtyp z2 = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
        
        numtyp recip = ucl_recip(r);
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip;
        fetch(psip,j,fp_tex);
        psip = ifp*rhojp + psip*rhoip + phip; 
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

__kernel void k_eam_fast(const __global numtyp4 *x_, 
                         const __global numtyp *fp_,
                         const __global int2 *type2rhor_z2r_in,
                         const __global numtyp4 *rhor_spline1, 
                         const __global numtyp4 *z2r_spline1,
                         const __global numtyp4 *z2r_spline2,
                         const __global int *dev_nbor, 
                         const __global int *dev_packed, 
                         __global acctyp4 *ans, 
                         __global acctyp *engv, 
                         const int eflag, const int vflag, const int inum, 
                         const int nbor_pitch, const numtyp cutforcesq, 
                         const numtyp rdr, const int nr, const int t_per_atom) {
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
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp ifp; fetch(ifp,i,fp_tex); //fp_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
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
        fetch4(coeff,index,rhor_sp1_tex);
        numtyp rhoip = (coeff.x*p + coeff.y)*p + coeff.z;
        
        mtype = jtype+iw;
        index = type2rhor_z2r[mtype].x*(nr+1)+m;
        fetch4(coeff,index,rhor_sp1_tex);
        numtyp rhojp = (coeff.x*p + coeff.y)*p + coeff.z;
        
        mtype = itype+jw;
        index = type2rhor_z2r[mtype].y*(nr+1)+m;
        fetch4(coeff,index,z2r_sp1_tex);
        numtyp z2p = (coeff.x*p + coeff.y)*p + coeff.z;
        fetch4(coeff,index,z2r_sp2_tex);
        numtyp z2 = ((coeff.x*p + coeff.y)*p + coeff.z)*p + coeff.w;
      
        numtyp recip = ucl_recip(r);
        numtyp phi = z2*recip;
        numtyp phip = z2p*recip - phi*recip;
        numtyp psip;
        fetch(psip,j,fp_tex);
        psip = ifp*rhojp + psip*rhoip + phip; 
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

