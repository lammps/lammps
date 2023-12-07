// **************************************************************************
//                                   edpd.cu
//                             -------------------
//                           Trung Dac Nguyen (U Chicago)
//
//  Device code for acceleration of the edpd pair style
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

#define store_drhoE(drhoI, deltaE, ii, inum, tid, t_per_atom, offset,        \
        drhoE)                                                               \
  if (t_per_atom>1) {                                                        \
    simdsync();                                                              \
    simd_reduce_add2(t_per_atom, red_acc, offset, tid, rhoEi, deltaE);       \
  }                                                                          \
  if (offset==0 && ii<inum) {                                                \
    drhoE[ii].x=drhoI;                                                       \
    drhoE[ii].y=deltaE;                                                      \
  }
#else
#define store_drhoE(drhoI, deltaE, ii, inum, tid, t_per_atom, offset, drhoE) \
  if (t_per_atom>1) {                                                        \
    simd_reduce_add2(t_per_atom,drhoI,deltaE);                               \
  }                                                                          \
  if (offset==0 && ii<inum) {                                                \
    drhoE[ii].x=drhoI;                                                       \
    drhoE[ii].y=deltaE;                                                      \
  }
#endif

__kernel void k_sph_lj(const __global numtyp4 *restrict x_,
                     const __global numtyp4 *restrict extra,
                     const __global numtyp4 *restrict coeff,
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
  acctyp Qi = (acctyp)0;

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp mass_itype = mass[itype];
    numtyp4 iv; fetch4(iv,i,vel_tex); //v_[i];
    int itag=iv.w;

    const numtyp4 Tcvi = extra[i];
    numtyp Ti = Tcvi.x;
    numtyp cvi = Tcvi.y;

    numtyp factor_dpd;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      factor_dpd = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      numtyp4 jv; fetch4(jv,j,vel_tex); //v_[j];
      int jtag=jv.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype]) {
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;
        numtyp vijeij = dot*rinv;

        const numtyp coeffx=coeff[mtype].x; // a0[itype][jtype]
        const numtyp coeffy=coeff[mtype].y; // gamma[itype][jtype]
        const numtyp coeffz=coeff[mtype].z; // cut[itype][jtype]

        const numtyp4 Tcvj = extra[j];
        numtyp Tj = Tcvj.x;
        numtyp cvj = Tcvj.y;

        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }

        numtyp randnum = (numtyp)0.0;
        saru(tag1, tag2, seed, timestep, randnum);

        numtyp T_ij=(numtyp)0.5*(Ti+Tj);
        numtyp4 T_pow;
        T_pow.x = T_ij - (numtyp)1.0;
        T_pow.y = T_pow.x*T_pow.x;
        T_pow.z = T_pow.x*T_pow.y;
        T_pow.w = T_pow.x*T_pow.z;

        numtyp coeff2x = coeff2[mtype].x; //power[itype][jtype]
        numtyp coeff2y = coeff2[mtype].y; //kappa[itype][jtype]
        numtyp coeff2z = coeff2[mtype].z; //powerT[itype][jtype]
        numtyp coeff2w = coeff2[mtype].w; //cutT[itype][jtype]
        numtyp power_d = coeff2x;
        if (power_flag) {
          numtyp factor = (numtyp)1.0;
          factor += sc[mtype].x*T_pow.x + sc[mtype].y*T_pow.y +
            sc[mtype].z*T_pow.z + sc[mtype].w*T_pow.w;
          power_d *= factor;
        }

        power_d = MAX((numtyp)0.01,power_d);
        numtyp wc = (numtyp)1.0 - r/coeffz; // cut[itype][jtype]
        wc = MAX((numtyp)0.0,MIN((numtyp)1.0,wc));
        numtyp wr = ucl_pow(wc, (numtyp)0.5*power_d);

        numtyp kboltz = (numtyp)1.0;
        numtyp GammaIJ = coeffy; // gamma[itype][jtype]
        numtyp SigmaIJ = (numtyp)4.0*GammaIJ*kboltz*Ti*Tj/(Ti+Tj);
        SigmaIJ = ucl_sqrt(SigmaIJ);

        numtyp force =  coeffx*T_ij*wc; // a0[itype][jtype]
        force -= GammaIJ *wr*wr *dot*rinv;
        force += SigmaIJ * wr *randnum * dtinvsqrt;
        force *= factor_dpd*rinv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        // heat transfer
        
        if (r < coeff2w) {  
          numtyp wrT = (numtyp)1.0 - r/coeff2w;
          wrT = MAX((numtyp)0.0,MIN((numtyp)1.0,wrT));
          wrT = ucl_pow(wrT, (numtyp)0.5*coeff2z); // powerT[itype][jtype]
          numtyp randnumT = (numtyp)0;
          saru(tag1, tag2, seed+tag1+tag2, timestep, randnumT); // randomT->gaussian();
          randnumT = MAX((numtyp)-5.0,MIN(randnum,(numtyp)5.0));

          numtyp kappaT = coeff2y; // kappa[itype][jtype]
          if (kappa_flag) {
            numtyp factor = (numtyp)1.0;
            factor += kc[mtype].x*T_pow.x + kc[mtype].y*T_pow.y +
              kc[mtype].z*T_pow.z + kc[mtype].w*T_pow.w;
            kappaT *= factor;
          }

          numtyp kij = cvi*cvj*kappaT * T_ij*T_ij;
          numtyp alphaij = ucl_sqrt((numtyp)2.0*kboltz*kij);

          numtyp dQc = kij * wrT*wrT * (Tj - Ti)/(Ti*Tj);
          numtyp dQd = wr*wr*( GammaIJ * vijeij*vijeij - SigmaIJ*SigmaIJ/mass_itype ) - SigmaIJ * wr *vijeij *randnum;
          dQd /= (cvi+cvj);
          numtyp dQr = alphaij * wrT * dtinvsqrt * randnumT;
          Qi += (dQc + dQd + dQr );
        }

        if (EVFLAG && eflag) {
          numtyp e = (numtyp)0.5*coeffx*T_ij*coeffz * wc*wc;
          energy+=factor_dpd*e;
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
  store_drhoE(Qi,ii,inum,tid,t_per_atom,offset,Q);
}

__kernel void k_sph_lj_fast(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict extra,
                          const __global numtyp2 *restrict coeff_in,
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
  const numtyp coeffy=coeff_in[ONETYPE].y;   // cutsq[itype][jtype]
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
  acctyp Qi = (acctyp)0;

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

    const numtyp4 Tcvi = extra[i];
    numtyp Ti = Tcvi.x;
    numtyp cvi = Tcvi.y;

    #ifndef ONETYPE
    numtyp factor_dpd;
    #endif
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
        numtyp r=ucl_sqrt(rsq);
        if (r < EPSILON) continue;

        numtyp rinv=ucl_recip(r);
        numtyp delvx = iv.x - jv.x;
        numtyp delvy = iv.y - jv.y;
        numtyp delvz = iv.z - jv.z;
        numtyp dot = delx*delvx + dely*delvy + delz*delvz;
        numtyp vijeij = dot*rinv;

        #ifndef ONETYPE
        const numtyp coeffx=coeff[mtype].x;   // a0[itype][jtype]
        const numtyp coeffy=coeff[mtype].y;   // gamma[itype][jtype]
        #endif

        const numtyp4 Tcvj = extra[j];
        numtyp Tj = Tcvj.x;
        numtyp cvj = Tcvj.y;

        unsigned int tag1=itag, tag2=jtag;
        if (tag1 > tag2) {
          tag1 = jtag; tag2 = itag;
        }
        numtyp randnum = (numtyp)0.0;
        saru(tag1, tag2, seed, timestep, randnum);

        numtyp T_ij=(numtyp)0.5*(Ti+Tj);
        numtyp4 T_pow;
        T_pow.x = T_ij - (numtyp)1.0;
        T_pow.y = T_pow.x*T_pow.x;
        T_pow.z = T_pow.x*T_pow.y;
        T_pow.w = T_pow.x*T_pow.z;

        numtyp power_d = coeff2x; // power[itype][jtype]
        if (power_flag) {
          numtyp factor = (numtyp)1.0;
          factor += scx*T_pow.x + scy*T_pow.y + scz*T_pow.z + scw*T_pow.w;
          power_d *= factor;
        }

        power_d = MAX((numtyp)0.01,power_d);
        numtyp wc = (numtyp)1.0 - r/coeffz; // cut[itype][jtype]
        wc = MAX((numtyp)0.0,MIN((numtyp)1.0,wc));
        numtyp wr = ucl_pow((numtyp)wc, (numtyp)0.5*power_d);

        numtyp kboltz = (numtyp)1.0;
        numtyp GammaIJ = coeffy; // gamma[itype][jtype]
        numtyp SigmaIJ = (numtyp)4.0*GammaIJ*kboltz*Ti*Tj/(Ti+Tj);
        SigmaIJ = ucl_sqrt(SigmaIJ);

        numtyp force =  coeffx*T_ij*wc; // a0[itype][jtype]
        force -= GammaIJ *wr*wr *dot*rinv;
        force += SigmaIJ* wr *randnum * dtinvsqrt;
        #ifndef ONETYPE
        force *= factor_dpd*rinv;
        #else
        force *= rinv;
        #endif

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        // heat transfer

        if (r < coeff2w) {  
          numtyp wrT = (numtyp)1.0 - r/coeff2w;
          wrT = MAX((numtyp)0.0,MIN((numtyp)1.0,wrT));
          wrT = ucl_pow(wrT, (numtyp)0.5*coeff2z); // powerT[itype][jtype]
          numtyp randnumT = (numtyp)0;
          saru(tag1, tag2, seed+tag1+tag2, timestep, randnumT); // randomT->gaussian();
          randnumT = MAX((numtyp)-5.0,MIN(randnum,(numtyp)5.0));

          numtyp kappaT = coeff2y; // kappa[itype][jtype]
          if (kappa_flag) {
            numtyp factor = (numtyp)1.0;
            factor += kcx*T_pow.x +  kcy*T_pow.y + kcz*T_pow.z + kcw*T_pow.w;
            kappaT *= factor;
          }
          
          numtyp kij = cvi*cvj*kappaT * T_ij*T_ij;
          numtyp alphaij = ucl_sqrt((numtyp)2.0*kboltz*kij);
         
          numtyp dQc = kij * wrT*wrT * (Tj - Ti )/(Ti*Tj);
          numtyp dQd = wr*wr*( GammaIJ * vijeij*vijeij - SigmaIJ*SigmaIJ/mass_itype ) - SigmaIJ * wr *vijeij *randnum;
          dQd /= (cvi+cvj);
          numtyp dQr = alphaij * wrT * dtinvsqrt * randnumT;
          Qi += (dQc + dQd + dQr );
        }

        if (EVFLAG && eflag) {
          numtyp e = (numtyp)0.5*coeffx*T_ij*coeffz * wc*wc;
          #ifndef ONETYPE
          energy+=factor_dpd*e;
          #else
          energy+=e;
          #endif
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
  store_drhoE(Qi,ii,inum,tid,t_per_atom,offset,Q);
}

