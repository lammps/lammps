// **************************************************************************
//                              sph_taitwater.cu
//                             -------------------
//                           Trung Dac Nguyen (U Chicago)
//
//  Device code for acceleration of the sph/taitwater pair style
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

#define store_drhoE(drhoEacc, ii, inum, tid, t_per_atom, offset, i, drhoE)  \
  if (t_per_atom>1) {                                                       \
    simdsync();                                                             \
    simd_reduce_add2(t_per_atom, red_acc, offset, tid,                      \
                     drhoEacc.x, drhoEacc.y);                               \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    drhoE[i]=drhoEacc.x;                                                    \
    drhoE[i+inum]=drhoEacc.y;                                               \
  }
#else
#define store_drhoE(drhoEacc, ii, inum, tid, t_per_atom, offset, i, drhoE)  \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      drhoEacc.x += shfl_down(drhoEacc.x, s, t_per_atom);                   \
      drhoEacc.y += shfl_down(drhoEacc.y, s, t_per_atom);                   \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    drhoE[i]=drhoEacc.x;                                                    \
    drhoE[i+inum]=drhoEacc.y;                                               \
  }
#endif

__kernel void k_sph_taitwater(const __global numtyp4 *restrict x_,
                              const __global numtyp4 *restrict extra,
                              const __global numtyp4 *restrict coeff,
                              const __global numtyp4 *restrict coeff2,
                              const int lj_types,
                              const __global numtyp *restrict sp_lj,
                              const __global int * dev_nbor,
                              const __global int * dev_packed,
                              __global acctyp3 *restrict ans,
                              __global acctyp *restrict engv,
                              __global acctyp *restrict drhoE,
                              const int eflag, const int vflag,
                              const int inum, const int nbor_pitch,
                              const __global numtyp4 *restrict v_,
                              const int dimension, const int t_per_atom) {
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }
  acctyp2 drhoEacc;
  drhoEacc.x = drhoEacc.y = (acctyp)0;

  if (ii<inum) {
    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp mass_itype = coeff2[itype].x;
    numtyp rho0_itype = coeff2[itype].y;
    numtyp soundspeed_itype = coeff2[itype].z;
    numtyp B_itype = coeff2[itype].w;
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;

    // compute pressure of atom i with Tait EOS
    numtyp tmp = rhoi / rho0_itype;
    numtyp fi = tmp * tmp * tmp;
    fi = B_itype * (fi * fi * tmp - (numtyp)1.0);
    fi /= (rhoi * rhoi);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<coeff[mtype].z) { // cutsq[itype][jtype]
        const numtyp coeffx=coeff[mtype].x;  // viscosity[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]

        numtyp mass_jtype = coeff2[jtype].x;
        numtyp rho0_jtype = coeff2[jtype].y;
        numtyp soundspeed_jtype = coeff2[jtype].z;
        numtyp B_jtype = coeff2[jtype].w;

        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;

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

        // compute pressure  of atom j with Tait EOS

        numtyp tmp = rhoj / rho0_jtype;
        numtyp fj = tmp * tmp * tmp;
        fj = B_jtype * (fj * fj * tmp - (numtyp)1.0);
        fj /= (rhoj * rhoj);

        // dot product of velocity delta and distance vector
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp delVdotDelR = delx*delvx + dely*delvy + delz*delvz;

        // artificial viscosity (Monaghan 1992)
        numtyp fvisc = (numtyp)0;
        if (delVdotDelR < (numtyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numtyp)0.01 * h * h);
          fvisc = -coeffx * (soundspeed_itype + soundspeed_jtype) * mu / (rhoi + rhoj);
        } 

        // total pair force & thermal energy increment
        numtyp force = -mass_itype * mass_jtype * (fi + fj + fvisc) * wfd;
        numtyp deltaE = (numtyp)-0.5 * force * delVdotDelR;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        // and change in density, drho[i]
        drhoEacc.x += mass_jtype * delVdotDelR * wfd;

        // change in thermal energy, desph[i]
        drhoEacc.y += deltaE;

        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }
    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
  store_drhoE(drhoEacc,ii,inum,tid,t_per_atom,offset,i,drhoE);
}

__kernel void k_sph_taitwater_fast(const __global numtyp4 *restrict x_,
                                   const __global numtyp4 *restrict extra,
                                   const __global numtyp4 *restrict coeff_in,
                                   const __global numtyp4 *restrict coeff2_in,
                                   const __global numtyp *restrict sp_lj_in,
                                   const __global int * dev_nbor,
                                   const __global int * dev_packed,
                                   __global acctyp3 *restrict ans,
                                   __global acctyp *restrict engv,
                                   __global acctyp *restrict drhoE,
                                   const int eflag, const int vflag,
                                   const int inum, const int nbor_pitch,
                                   const __global numtyp4 *restrict v_,
                                   const int dimension, const int t_per_atom) {
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  #ifndef ONETYPE
  __local numtyp4 coeff[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 coeff2[MAX_SHARED_TYPES];
  if (tid<MAX_SHARED_TYPES) {
    coeff2[tid] = coeff2_in[tid];
  }
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff[tid]=coeff_in[tid];
  }
  __syncthreads();
  #else
  const numtyp coeffx=coeff_in[ONETYPE].x;   // viscosity[itype][jtype]
  const numtyp coeffy=coeff_in[ONETYPE].y;   // cut[itype][jtype]
  const numtyp cutsq_p=coeff_in[ONETYPE].z;  // cutsq[itype][jtype]
  #endif

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }
  acctyp2 drhoEacc;
  drhoEacc.x = drhoEacc.y = (acctyp)0;

  if (ii<inum) {
    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    numtyp mass_itype = coeff2[iw].x;
    numtyp rho0_itype = coeff2[iw].y;
    numtyp soundspeed_itype = coeff2[iw].z;
    numtyp B_itype = coeff2[iw].w;
    #ifndef ONETYPE
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    #endif
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;

    // compute pressure of atom i with Tait EOS
    numtyp tmp = rhoi / rho0_itype;
    numtyp fi = tmp * tmp * tmp;
    fi = B_itype * (fi * fi * tmp - (numtyp)1.0);
    fi /= (rhoi * rhoi);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      #ifndef ONETYPE
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      const numtyp cutsq_p=coeff[mtype].z;
      #endif
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq_p) {
        #ifndef ONETYPE
        const numtyp coeffx=coeff[mtype].x;  // viscosity[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]
        #endif

        numtyp mass_jtype = coeff2[jtype].x;
        numtyp rho0_jtype = coeff2[jtype].y;
        numtyp soundspeed_jtype = coeff2[jtype].z;
        numtyp B_jtype = coeff2[jtype].w;

        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;

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

        // compute pressure  of atom j with Tait EOS
        numtyp tmp = rhoj / rho0_jtype;
        numtyp fj = tmp * tmp * tmp;
        fj = B_jtype * (fj * fj * tmp - (numtyp)1.0);
        fj /= (rhoj * rhoj);

        // dot product of velocity delta and distance vector
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp delVdotDelR = delx*delvx + dely*delvy + delz*delvz;

        // artificial viscosity (Monaghan 1992)
        numtyp fvisc = (numtyp)0;
        if (delVdotDelR < (numtyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numtyp)0.01 * h * h);
          fvisc = -coeffx * (soundspeed_itype + soundspeed_jtype) * mu / (rhoi + rhoj);
        }

        // total pair force & thermal energy increment
        numtyp force = -mass_itype * mass_jtype * (fi + fj + fvisc) * wfd;
        numtyp deltaE = (numtyp)-0.5 * force * delVdotDelR;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        // and change in density, drho[i]
        drhoEacc.x += mass_jtype * delVdotDelR * wfd;

        // change in thermal energy, desph[i]
        drhoEacc.y += deltaE;

        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }
    } // for nbor
  } // if ii

  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag, ans,engv);
  store_drhoE(drhoEacc,ii,inum,tid,t_per_atom,offset,i,drhoE);
}

