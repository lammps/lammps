// **************************************************************************
//                                 sph_lj.cu
//                             -------------------
//                           Trung Dac Nguyen (U Chicago)
//
//  Device code for acceleration of the sph/lj pair style
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
    drhoE[ii]=drhoEacc;                                                     \
  }
#endif

/* ------------------------------------------------------------------------ */
/* Lennard-Jones EOS,
   Francis H. Ree
   "Analytic representation of thermodynamic data for the Lennard‐Jones fluid",
   Journal of Chemical Physics 73 pp. 5401-5403 (1980)
   return p = pc[0], c = pc[1]
*/

ucl_inline void LJEOS2(const numtyp rho, const numtyp e, const numtyp cv, numtyp pc[2])
{
  numtyp T = e/cv;
  numtyp beta = ucl_recip(T); // (numtyp)1.0 / T;
  numtyp beta_sqrt = ucl_sqrt(beta);
  numtyp x = rho * ucl_sqrt(beta_sqrt);

  numtyp xsq = x * x;
  numtyp xpow3 = xsq * x;
  numtyp xpow4 = xsq * xsq;

  /* differential of Helmholtz free energy w.r.t. x */
  numtyp diff_A_NkT = (numtyp)3.629 + (numtyp)7.264*x -
    beta*((numtyp)3.492 - (numtyp)18.698*x + (numtyp)35.505*xsq - (numtyp)31.816*xpow3 +
    (numtyp)11.195*xpow4) - beta_sqrt*((numtyp)5.369 + (numtyp)13.16*x +
    (numtyp)18.525*xsq - (numtyp)17.076*xpow3 + (numtyp)9.32*xpow4) +
    (numtyp)10.4925*xsq + (numtyp)11.46*xpow3 + (numtyp)2.176*xpow4*xpow4*x;

  /* differential of Helmholtz free energy w.r.t. x^2 */
  numtyp d2A_dx2 = (numtyp)7.264 + (numtyp)20.985*x +
     beta*((numtyp)18.698 - (numtyp)71.01*x + (numtyp)95.448*xsq - (numtyp)44.78*xpow3) -
     beta_sqrt*((numtyp)13.16 + (numtyp)37.05*x - (numtyp)51.228*xsq + (numtyp)37.28*xpow3) +
     (numtyp)34.38*xsq + (numtyp)19.584*xpow4*xpow4;

  // p = rho k T * (1 + rho * d(A/(NkT))/drho)
  // dx/drho = rho/x
  pc[0] = rho * T * ((numtyp)1.0 + diff_A_NkT * x); // pressure
  numtyp csq = T * ((numtyp)1.0 + (numtyp)2.0 * diff_A_NkT * x + d2A_dx2 * x * x); // soundspeed squared
  if (csq > (numtyp)0.0) {
    pc[1] = ucl_sqrt(csq); // soundspeed
  } else {
    pc[1] = (numtyp)0.0;
  }
}


__kernel void k_sph_lj(const __global numtyp4 *restrict x_,
                       const __global numtyp4 *restrict extra,
                       const __global numtyp4 *restrict coeff,
                       const __global numtyp *restrict mass,
                       const int lj_types,
                       const __global numtyp *restrict sp_lj,
                       const __global int * dev_nbor,
                       const __global int * dev_packed,
                       __global acctyp3 *restrict ans,
                       __global acctyp *restrict engv,
                       __global acctyp2 *restrict drhoE,
                       const int eflag, const int vflag,
                       const int inum, const int nbor_pitch,
                       const __global numtyp4 *restrict v_,
                       const int dimension, const int t_per_atom) {
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
    numtyp mass_itype = mass[itype];
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp esphi = extrai.y;
    numtyp cvi = extrai.z;

    // compute pressure of particle i with LJ EOS
    numtyp fci[2];
    LJEOS2(rhoi, esphi, cvi, fci);
    numtyp fi = fci[0];
    numtyp ci = fci[1];
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
        numtyp mass_jtype = mass[jtype];
        const numtyp coeffx=coeff[mtype].x;  // viscosity[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]

        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;
        numtyp cvj = extraj.z;

        numtyp h = coeffy; // cut[itype][jtype]
        numtyp ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;
        numtyp ihcub = ihsq * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = (numtyp)-19.098593171027440292 * wfd * wfd * ihsq * ihsq * ihsq;
        }

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
        if (delVdotDelR < (numtyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numtyp)0.01 * h * h);
          fvisc = -coeffx * (ci + cj) * mu / (rhoi + rhoj); // viscosity[itype][jtype]
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

__kernel void k_sph_lj_fast(const __global numtyp4 *restrict x_,
                            const __global numtyp4 *restrict extra,
                            const __global numtyp4 *restrict coeff_in,
                            const __global numtyp *restrict mass,
                            const __global numtyp *restrict sp_lj_in,
                            const __global int * dev_nbor,
                            const __global int * dev_packed,
                            __global acctyp3 *restrict ans,
                            __global acctyp *restrict engv,
                            __global acctyp2 *restrict drhoE,
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

    const numtyp4 extrai = extra[i];
    numtyp rhoi = extrai.x;
    numtyp esphi = extrai.y;
    numtyp cvi = extrai.z;

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
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype = jx.w;
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      const numtyp cutsq_p=coeff[mtype].z; // cutsq[itype][jtype];
      #endif
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq_p) {
        numtyp mass_jtype = mass[jtype];
        #ifndef ONETYPE
        const numtyp coeffx=coeff[mtype].x;  // viscosity[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;  // cut[itype][jtype]
        #endif
        const numtyp4 extraj = extra[j];
        numtyp rhoj = extraj.x;
        numtyp esphj = extraj.y;
        numtyp cvj = extraj.z;

        numtyp h = coeffy; // cut[itype][jtype]
        numtyp ih = ucl_recip(h); // (numtyp)1.0 / h;
        numtyp ihsq = ih * ih;
        numtyp ihcub = ihsq * ih;

        numtyp wfd = h - ucl_sqrt(rsq);
        if (dimension == 3) {
          // Lucy Kernel, 3d
          wfd = (numtyp)-25.066903536973515383 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = (numtyp)-19.098593171027440292 * wfd * wfd * ihsq * ihsq * ihsq;
        }

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
        if (delVdotDelR < (numtyp)0) {
          numtyp mu = h * delVdotDelR / (rsq + (numtyp)0.01 * h * h);
          fvisc = -coeffx * (ci + cj) * mu / (rhoi + rhoj); // viscosity[itype][jtype]
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

