// **************************************************************************
//                             sph_heatconduction.cu
//                             ---------------------
//                           Trung Dac Nguyen (U Chicago)
//
//  Device code for acceleration of the sph/heatconduction pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : September 2023
//    email                : ndactrung@gmail.com
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( vel_tex,float4);
#else
_texture_2d( pos_tex,int4);
_texture_2d( vel_tex,int4);
#endif
#else
#define pos_tex x_
#define vel_tex v_
#endif

#if (SHUFFLE_AVAIL == 0)

#define store_dE(dEacc, ii, inum, tid, t_per_atom, offset, dE)              \
  if (t_per_atom>1) {                                                       \
    simdsync();                                                             \
    simd_reduce_add1(t_per_atom, red_acc, offset, tid, dEacc);              \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    dE[ii]=dEacc;                                                           \
  }
#else
#define store_drhoE(dEacc, ii, inum, tid, t_per_atom, offset, dE)           \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      dEacc += shfl_down(dEacc, s, t_per_atom);                             \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    dE[ii]=dEacc;                                                           \
  }
#endif

/* ------------------------------------------------------------------------ */

__kernel void k_sph_heatconduction(const __global numtyp4 *restrict x_,
                       const __global numtyp4 *restrict extra,
                       const __global numtyp4 *restrict coeff,
                       const __global numtyp *restrict mass,
                       const int lj_types,
                       const __global numtyp *restrict sp_lj,
                       const __global int * dev_nbor,
                       const __global int * dev_packed,
                       __global acctyp3 *restrict ans,
                       __global acctyp *restrict engv,
                       __global acctyp *restrict dE,
                       const int eflag, const int vflag,
                       const int inum, const int nbor_pitch,
                       const __global numtyp4 *restrict v_,
                       const int dimension, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
#if (SHUFFLE_AVAIL == 0)
  local_allocate_store_pair();
#endif

  acctyp dEacc = (acctyp)0;

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp mass_itype = mass[itype];

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp esphi = extrai.y;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<coeff[mtype].z) { // cutsq[itype][jtype]
        numtyp mass_jtype = mass[jtype];
        const numtyp coeffx=coeff[mtype].x;  // alpha[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]

        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;

        numtyp h = coeffy; // cut[itype][jtype]
        numtyp ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = (numtyp)-19.098593171027440292 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        // total thermal energy increment
        numtyp D = coeffx; // alpha[itype][jtype]  diffusion coefficient
        numtyp deltaE = (numtyp)2.0 * mass_itype * mass_jtype / (mass_itype + mass_jtype);
        deltaE *= (rhoi + rhoj) / (rhoi * rhoj);
        deltaE *= D * (esphi - esphj) * wfd;

        // change in thermal energy, desph[i]
        dEacc += deltaE;

      }
    } // for nbor
  } // if ii

  store_drhoE(dEacc,ii,inum,tid,t_per_atom,offset,dE);
}

__kernel void k_sph_heatconduction_fast(const __global numtyp4 *restrict x_,
                            const __global numtyp4 *restrict extra,
                            const __global numtyp4 *restrict coeff_in,
                            const __global numtyp *restrict mass,
                            const __global numtyp *restrict sp_lj_in,
                            const __global int * dev_nbor,
                            const __global int * dev_packed,
                            __global acctyp3 *restrict ans,
                            __global acctyp *restrict engv,
                            __global acctyp *restrict dE,
                            const int eflag, const int vflag,
                            const int inum, const int nbor_pitch,
                            const __global numtyp4 *restrict v_,
                            const int dimension, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  #ifndef ONETYPE
  __local numtyp4 coeff[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff[tid]=coeff_in[tid];
  }
  __syncthreads();
  #else
  const numtyp coeffx=coeff_in[ONETYPE].x;   // alpha[itype][jtype]
  const numtyp coeffy=coeff_in[ONETYPE].y;   // cut[itype][jtype]
  const numtyp cutsq_p=coeff_in[ONETYPE].z;  // cutsq[itype][jtype]
  #endif

  int n_stride;
#if (SHUFFLE_AVAIL == 0)
  local_allocate_store_pair();
#endif

  acctyp dEacc = (acctyp)0;

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    numtyp mass_itype = mass[iw];
    #ifndef ONETYPE
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    #endif

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp esphi = extrai.y;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      #ifndef ONETYPE
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype = jx.w;
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      const numtyp cutsq_p=coeff[mtype].z;
      #endif

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq_p) {
        numtyp mass_jtype = mass[jtype];
        #ifndef ONETYPE
        const numtyp coeffx=coeff[mtype].x;  // alpha[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]
        #endif
        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;

        numtyp h = coeffy; // cut[itype][jtype]
        numtyp ih = ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = (numtyp)-19.098593171027440292 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        // total thermal energy increment
        numtyp D = coeffx; // alpha[itype][jtype]  diffusion coefficient
        numtyp deltaE = (numtyp)2.0 * mass_itype * mass_jtype / (mass_itype + mass_jtype);
        deltaE *= (rhoi + rhoj) / (rhoi * rhoj);
        deltaE *= D * (esphi - esphj) * wfd;

        // change in thermal energy, desph[i]
        dEacc += deltaE;

      }
    } // for nbor
  } // if ii

  store_drhoE(dEacc,ii,inum,tid,t_per_atom,offset,dE);
}

