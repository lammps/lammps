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

#define store_drhoE(drhoEacc, ii, inum, tid, t_per_atom, offset, drhoE)      \
  if (t_per_atom>1) {                                                        \
    simdsync();                                                              \
    simd_reduce_add2(t_per_atom, red_acc, offset, tid,                       \
                     drhoEacc.x, drhoEacc.y);                                \
  }                                                                          \
  if (offset==0 && ii<inum) {                                                \
    drhoE[ii]=drhoEacc;                                                      \
  }
#else
#define store_drhoE(drhoEacc, ii, inum, tid, t_per_atom, offset, drhoE)     \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      drhoEacc.x += shfl_down(drhoEacc.x, s, t_per_atom);                   \
      drhoEacc.y += shfl_down(drhoEacc.y, s, t_per_atom);                   \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    drhoE[ii].x=drhoEacc;                                                   \
  }
#endif

__kernel void k_sph_taitwater(const __global numtyp4 *restrict x_,
                       const __global numtyp4 *restrict extra,
                       const __global numtyp4 *restrict coeff,
                       const __global numtyp2 *restrict rho0sspeed,
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
                       const int t_per_atom) {
  int tid, ii, offset;
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
  drhoEacc.x = drhoEacc.x = (acctyp)0;

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp massi= extrai.y;

    // compute pressure of atom i with Tait EOS
    numtyp tmp = rho[i] / rho0[itype];
    numtyp fi = tmp * tmp * tmp;
    fi = B[itype] * (fi * fi * tmp - 1.0);
    fi /= (rhoi * rhoi);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
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

        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;
        numtyp cvj = extraj.z;
        numtyp massj= extraj.w;

        numtyp h = coeffy; // cut[itype][jtype]
        ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;
        numtyp ihcub = ihsq * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        // domain->dimension == 3 Lucy Kernel, 3d
        wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
       
        // Lucy Kernel, 2d
        //wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        
        // compute pressure  of atom j with Tait EOS
        numtyp tmp = rho[j] / rho0[jtype];
        numtyp fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0);
        fj /= (rhoj * rhoj);

        // dot product of velocity delta and distance vector
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp delVdotDelR = delx*delvx + dely*delvy + delz*delvz;

        // artificial viscosity (Monaghan 1992)
        numtyp fvisc = (numtyp)0;
        if (delVdotDelR < (numyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numyp)0.01 * h * h);
          //fvisc = -coeffx * (ci + cj) * mu / (rhoi + rhoj); // viscosity[itype][jtype]

          fvisc = -coeffx * (soundspeed[itype]
              + soundspeed[jtype]) * mu / (rhoi + rhoj);
        }

        // total pair force & thermal energy increment
        numtyp force = -massi * massj * (fi + fj + fvisc) * wfd;
        numtyp deltaE = (numtyp)-0.5 * force * delVdotDelR;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        // and change in density, drho[i]
        drhoEacc.x += massj * delVdotDelR * wfd;

        // change in thermal energy, desph[i]
        drhoEacc.y += deltaE;

        if (EVFLAG && eflag) {
          numtyp e = (numtyp)0;
          energy+=e;
        }
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
  store_drhoE(drhoEacc,ii,inum,tid,t_per_atom,offset,drhoE);
}

__kernel void k_sph_taitwater_fast(const __global numtyp4 *restrict x_,
                            const __global numtyp4 *restrict extra,
                            const __global numtyp4 *restrict coeff_in,
                            const __global numtyp2 *restrict rho0sspeed_in,
                            const __global numtyp *restrict sp_lj_in,
                            const __global int * dev_nbor,
                            const __global int * dev_packed,
                            __global acctyp3 *restrict ans,
                            __global acctyp *restrict engv,
                            __global acctyp *restrict drhoE,
                            const int eflag, const int vflag,
                            const int inum, const int nbor_pitch,
                            const __global numtyp4 *restrict v_,
                            const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  #ifndef ONETYPE
  __local numtyp4 coeff[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4) {
    sp_lj[tid]=sp_lj_in[tid];
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
  drhoEacc.x = drhoEacc.x = (acctyp)0;

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
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp esphi = extrai.y;
    numtyp cvi = extrai.z;
    numtyp massi= extrai.w;

    // compute pressure of particle i with LJ EOS
    numtyp fci[2];
    LJEOS2(rhoi, esphi, cvi, fci);
    numtyp fi = fci[0];
    numtyp ci = fci[1];
    fi /= (rhoi * rhoi);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      #ifndef ONETYPE
      factor_dpd = sp_lj[sbmask(j)];
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      const numtyp cutsq_p=cutsq[mtype];
      #endif
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];
      int jtag=jv.w;

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
        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;
        numtyp cvj = extraj.z;
        numtyp massj= extraj.w;

        numtyp h = coeffy; // cut[itype][jtype]
        ih = ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;
        numtyp ihcub = ihsq * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        // domain->dimension == 3 Lucy Kernel, 3d
        wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
       
        // Lucy Kernel, 2d
        //wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        
        // function call to LJ EOS
        numtyp fcj[2];
        LJEOS2(rhoj, esphj, cvj, fcj);
        numtyp fj = fcj[0];
        numtyp cj = fcj[1];
        fj /= (rhoj * rhoj);

        // apply long-range correction to model a LJ fluid with cutoff
        // this implies that the modelled LJ fluid has cutoff == SPH cutoff
        numtyp lrc = (numtyp)-11.1701 * (ihcub * ihcub * ihcub - (numtyp)1.5 * ihcub);
        fi += lrc;
        fj += lrc;

        // dot product of velocity delta and distance vector
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp delVdotDelR = delx*delvx + dely*delvy + delz*delvz;

        // artificial viscosity (Monaghan 1992)
        numtyp fvisc = (numtyp)0;
        if (delVdotDelR < (numyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numyp)0.01 * h * h);
          fvisc = -coeffx * (ci + cj) * mu / (rhoi + rhoj); // viscosity[itype][jtype]
        }

        // total pair force & thermal energy increment
        numtyp force = -massi * massj * (fi + fj + fvisc) * wfd;
        numtyp deltaE = (numtyp)-0.5 * force * delVdotDelR;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          numtyp e = (numtyp)0;          
          energy+=e;
        }
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
  store_drhoE(drhoEacc,ii,inum,tid,t_per_atom,offset,drhoE);
}

